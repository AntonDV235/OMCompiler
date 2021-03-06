/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */

#include <string.h>
#include <setjmp.h>

#include "openmodelica.h"
#include "openmodelica_func.h"
#include "simulation_data.h"

#include "util/omc_error.h"
#include "util/memory_pool.h"

#include "simulation/options.h"
#include "simulation/simulation_runtime.h"
#include "simulation/results/simulation_result.h"
#include "simulation/solver/solver_main.h"
#include "simulation/solver/model_help.h"
#include "simulation/solver/external_input.h"
#include "simulation/solver/epsilon.h"

#include "simulation/solver/dassl.h"
#include "meta/meta_modelica.h"

#ifdef __cplusplus
extern "C" {
#endif

static const char *dasslJacobianMethodStr[DASSL_JAC_MAX] = {"unknown",
                                                "coloredNumerical",
                                                "coloredSymbolical",
                                                "internalNumerical",
                                                "numerical",
                                                "symbolical"
                                                };

static const char *dasslJacobianMethodDescStr[DASSL_JAC_MAX] = {"unknown",
                                                       "colored numerical jacobian - default.",
                                                       "colored symbolic jacobian - needs omc compiler flags +generateSymbolicJacobian or +generateSymbolicLinearization.",
                                                       "internal numerical jacobian.",
                                                       "numerical jacobian.",
                                                       "symbolic jacobian - needs omc compiler flags +generateSymbolicJacobian or +generateSymbolicLinearization."
                                                      };

/* experimental flag for SKF TLM Master Solver Interface
 *  - it's used with -noEquidistantTimeGrid flag.
 *  - it's set to 1 if the continuous system is evaluated
 *    when dassl finished a step, otherwise it's 0.
 */
int RHSFinalFlag;

/* provides a dummy Jacobian to be used with DASSL */
static int
dummy_Jacobian(double *t, double *y, double *yprime, double *deltaD,
    double *delta, double *cj, double *h, double *wt, double *rpar, int* ipar) {
  return 0;
}
static int
dummy_zeroCrossing(int *neqm, double *t, double *y, double *yp,
                   int *ng, double *gout, double *rpar, int* ipar) {
  return 0;
}

static int JacobianSymbolic(double *t, double *y, double *yprime,  double *deltaD, double *pd, double *cj, double *h, double *wt,
    double *rpar, int* ipar);
static int JacobianSymbolicColored(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
    double *rpar, int* ipar);
static int JacobianOwnNum(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
    double *rpar, int* ipar);
static int JacobianOwnNumColored(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
    double *rpar, int* ipar);

void  DDASKR(
    int (*res) (double *t, double *y, double *yprime, double* cj, double *delta, int *ires, double *rpar, int* ipar),
    int *neq,
    double *t,
    double *y,
    double *yprime,
    double *tout,
    int *info,
    double *rtol,
    double *atol,
    int *idid,
    double *rwork,
    int *lrw,
    int *iwork,
    int *liw,
    double *rpar,
    int *ipar,
    int (*jac) (double *t, double *y, double *yprime, double *deltaD, double *delta, double *cj, double *h, double *wt, double *rpar, int* ipar),
    int (*psol) (int *neq, double *t, double *y, double *yprime, double *savr, double *pwk, double *cj, double *wt, double *wp, int *iwp, double *b, double eplin, int* ires, double *rpar, int* ipar),
    int (*g) (int *neqm, double *t, double *y, double *yp, int *ng, double *gout, double *rpar, int* ipar),
    int *ng,
    int *jroot
);

static int
dummy_precondition(int *neq, double *t, double *y, double *yprime, double *savr, double *pwk, double *cj, double *wt, double *wp, int *iwp, double *b, double eplin, int* ires, double *rpar, int* ipar){
    return 0;
}

static int continue_DASSL(int* idid, double* tolarence);

/* function for calculating state values on residual form */
static int functionODE_residual(double *t, double *x, double *xprime, double *cj, double *delta, int *ires, double *rpar, int* ipar);
/* function for calculating zeroCrossings */
static int function_ZeroCrossingsDASSL(int *neqm, double *t, double *y, double *yp,
        int *ng, double *gout, double *rpar, int* ipar);

int dassl_initial(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo, DASSL_DATA *dasslData)
{
  TRACE_PUSH
  /* work arrays for DASSL */
  unsigned int i;
  SIMULATION_DATA tmpSimData = {0};

  RHSFinalFlag = 0;

  dasslData->liw = 40 + data->modelData->nStates;
  dasslData->lrw = 60 + ((maxOrder + 4) * data->modelData->nStates) + (data->modelData->nStates * data->modelData->nStates)  + (3*data->modelData->nZeroCrossings);
  dasslData->rwork = (double*) calloc(dasslData->lrw, sizeof(double));
  assertStreamPrint(threadData, 0 != dasslData->rwork,"out of memory");
  dasslData->iwork = (int*)  calloc(dasslData->liw, sizeof(int));
  assertStreamPrint(threadData, 0 != dasslData->iwork,"out of memory");
  dasslData->ng = (int) data->modelData->nZeroCrossings;
  dasslData->jroot = (int*)  calloc(data->modelData->nZeroCrossings, sizeof(int));
  dasslData->rpar = (double**) malloc(3*sizeof(double*));
  dasslData->ipar = (int*) malloc(sizeof(int));
  dasslData->ipar[0] = ACTIVE_STREAM(LOG_JAC);
  assertStreamPrint(threadData, 0 != dasslData->ipar,"out of memory");
  dasslData->atol = (double*) malloc(data->modelData->nStates*sizeof(double));
  dasslData->rtol = (double*) malloc(data->modelData->nStates*sizeof(double));
  dasslData->info = (int*) calloc(infoLength, sizeof(int));
  assertStreamPrint(threadData, 0 != dasslData->info,"out of memory");
  dasslData->dasslStatistics = (unsigned int*) calloc(numStatistics, sizeof(unsigned int));
  assertStreamPrint(threadData, 0 != dasslData->dasslStatistics,"out of memory");
  dasslData->dasslStatisticsTmp = (unsigned int*) calloc(numStatistics, sizeof(unsigned int));
  assertStreamPrint(threadData, 0 != dasslData->dasslStatisticsTmp,"out of memory");

  dasslData->idid = 0;

  dasslData->sqrteps = sqrt(DBL_EPSILON);
  dasslData->ysave = (double*) malloc(data->modelData->nStates*sizeof(double));
  dasslData->delta_hh = (double*) malloc(data->modelData->nStates*sizeof(double));
  dasslData->newdelta = (double*) malloc(data->modelData->nStates*sizeof(double));
  dasslData->stateDer = (double*) malloc(data->modelData->nStates*sizeof(double));

  data->simulationInfo->currentContext = CONTEXT_ALGEBRAIC;

  /* ### start configuration of dassl ### */
  infoStreamPrint(LOG_SOLVER, 1, "Configuration of the dassl code:");



  printf("---------------###---------------\n\n");
  printf("Some stats:\n\n");
  printf("(A) Firstly, we consider some info of the object DATA *data:\n");
  printf("(1) SIMULATION_INFO (*simulationInfo):\n");
  printf("\tstartTime: %f\n", data->simulationInfo->startTime);
  printf("\tstopTime: %f\n", data->simulationInfo->stopTime);
  printf("\tnumSteps: %d\n", data->simulationInfo->numSteps);
  printf("\ttolerance: %f\n", data->simulationInfo->tolerance);
  printf("\ttolerance: %f\n", (*(*data).simulationInfo).tolerance);
  printf("\t*solverMethod: %c\n", *(data->simulationInfo->solverMethod));
  printf("\t*outputFormat: %c\n", *(data->simulationInfo->outputFormat));
  printf("\t*variableFilter: %c\n", *(data->simulationInfo->variableFilter));
  printf("\tlsMethod: %d\n", data->simulationInfo->lsMethod);
  printf("\tmixedMethod: %d\n", data->simulationInfo->mixedMethod);
  printf("\tnlsMethod: %d\n", data->simulationInfo->nlsMethod);
  printf("\tnewtonStrategy: %d\n", data->simulationInfo->newtonStrategy);
  printf("\tnlsCsvInfomation: %d\n", data->simulationInfo->nlsCsvInfomation);
  printf("\tcurrentContext: %d\n", data->simulationInfo->currentContext);
  printf("\tcurrentContextOld: %d\n", data->simulationInfo->currentContextOld);
  printf("\tjacobianEvals: %d\n", data->simulationInfo->jacobianEvals);
  printf("\tlambda: %d\n", data->simulationInfo->lambda);
  printf("\tnextSampleEvent: %f\n", data->simulationInfo->nextSampleEvent);
  printf("\tnextSampleTimes[0]: %f\n", data->simulationInfo->nextSampleTimes[0]);
  printf("\tnextSampleTimes[1]: %f\n", data->simulationInfo->nextSampleTimes[1]);
  printf("\tnextSampleTimes[2]: %f\n", data->simulationInfo->nextSampleTimes[2]);


  printf("\n\tThere are even more although I believe they will not be useful for now.\n\n");

  printf("(2) SIMULATION_DATA (**localData):\n");
  printf("\ttimeValue: %f\n", data->localData[0]->timeValue);
  printf("The following are the initial values.\n");
  printf("\trealVars[0]: %f\n", data->localData[0]->realVars[0]);
  printf("\trealVars[1]: %f\n", data->localData[0]->realVars[1]);
  printf("\tintegerVars: %d\n", data->localData[0]->integerVars[1]);
  printf("\t*booleanVars: %s\n", *(data->localData[0]->booleanVars));
  printf("\t*stringVars: %s\n", *(data->localData[0]->stringVars));
  //printf("\t*inlineVars: %f\n", *(data->localData[0]->inlineVars)); //Errors when attempting to print this one. Not sure why though.
  printf("\n\tI am not sure this garbage struct is necessary. One \'timeValue\' is used. I am sure we can use a local variable or a similar route.\n\n");

  printf("(3) MODEL_DATA (*modelData):\n");
  printf("\trealVarsData[0].attribute.nominal: %f\n", data->modelData->realVarsData[0].attribute.nominal);
  printf("\trealVarsData[1].attribute.nominal: %f\n", data->modelData->realVarsData[1].attribute.nominal);
  printf("\trealVarsData[0].attribute.useStart: %s\n", data->modelData->realVarsData[0].attribute.useStart ? "true" : "false");
  printf("\trealVarsData[1].attribute.useStart: %s\n", data->modelData->realVarsData[1].attribute.useStart ? "true" : "false");
  printf("\tnStates: %d\n", data->modelData->nStates);
  printf("\tnVariablesReal: %d\n", data->modelData->nVariablesReal);
  printf("\tnDiscreteReal: %d\n", data->modelData->nDiscreteReal);
  printf("\tnVariablesInteger: %d\n", data->modelData->nVariablesInteger);
  printf("\tnVariablesBoolean: %d\n", data->modelData->nVariablesBoolean);
  printf("\tnVariablesString: %d\n", data->modelData->nVariablesString);
  printf("\tnParametersReal: %d\n", data->modelData->nParametersReal);
  printf("\tnVariablesReal: %d\n", data->modelData->nVariablesReal);
  printf("\tnParametersInteger: %d\n", data->modelData->nParametersInteger);
  printf("\tnParametersBoolean: %d\n", data->modelData->nParametersBoolean);
  printf("\tnParametersString: %d\n", data->modelData->nParametersString);
  printf("\tnInputVars: %d\n", data->modelData->nInputVars);
  printf("\tnOutputVars: %d\n", data->modelData->nOutputVars);
  printf("\tnZeroCrossings: %d\n", data->modelData->nZeroCrossings);
  printf("\tnMathEvents: %d\n", data->modelData->nMathEvents);
  printf("\tnDelayExpressions: %d\n", data->modelData->nDelayExpressions);
  printf("\tnExtObjs: %d\n", data->modelData->nExtObjs);
  printf("\tnMixedSystems: %d\n", data->modelData->nMixedSystems);
  printf("\tnLinearSystems: %d\n", data->modelData->nLinearSystems);
  printf("\tnNonLinearSystems: %d\n", data->modelData->nNonLinearSystems);
  printf("\tnStateSets: %d\n", data->modelData->nStateSets);
  printf("\tnInlineVars: %d\n", data->modelData->nInlineVars);
  printf("\tnOptimizeConstraints: %d\n", data->modelData->nOptimizeConstraints);
  printf("\tnOptimizeFinalConstraints: %d\n\n", data->modelData->nOptimizeFinalConstraints);

  printf("(4) RINGBUFFER* (simulationData) -- The first comment in the .c file is that this does not work. This struct is never used in the qss solver.\n\n");
  //printf("\titemSize: %d\n", data->simulationData->itemSize);

  printf("(5) OpenModelicaGeneratedFunctionCallbacks (*callback). Just a number of functions updating the results. Not even sure where all these functions are defined. We may want to investigate where these functions are defined if we want to have a full understanding of the how OpenModelica enables it's solvers.\n\n");

  printf("(B) threadData_t threadData. Next we have to consider the object threadData_t. However, I cannot find the definition of this object. Must investigate and get all attributes.\n\n");

  printf("(C) SOLVER_INFO solverInfo:\n");
  printf("\tcurrentTime: %f\n", solverInfo->currentTime);
  printf("\tcurrentStepSize: %f\n", solverInfo->currentStepSize);
  printf("\tlaststep: %f\n", solverInfo->laststep);
  printf("\tsolverMethod: %d\n", solverInfo->solverMethod);
  printf("\tsolverStepSize: %f\n", solverInfo->solverStepSize);
  printf("\tsolverRootFinding: %s\n", solverInfo->solverRootFinding ? "true" : "false");
  printf("\tsolverNoEquidistantGrid: %s\n", solverInfo->solverNoEquidistantGrid ? "true" : "false");
  printf("\tlastdesiredStep: %f\n", solverInfo->lastdesiredStep);
  printf("\tstateEvents: %d\n", solverInfo->stateEvents);
  printf("\tdidEventStep: %d\n", solverInfo->didEventStep);
//  printf("\teventLst.itemSize %d\n", solverInfo->eventLst->itemSize);
//  printf("\teventLst.length %d\n", solverInfo->eventLst.length);
  printf("\tsampleEvents: %d\n", solverInfo->sampleEvents);
  printf("\tintegratorSteps: %d\n", solverInfo->integratorSteps);


  printf("---------------###---------------\n\n");


  /* set nominal values of the states for absolute tolerances */
  dasslData->info[1] = 1;
  infoStreamPrint(LOG_SOLVER, 1, "The relative tolerance is %g. Following absolute tolerances are used for the states: ", data->simulationInfo->tolerance);
  for(i=0; i<data->modelData->nStates; ++i)
  {
    dasslData->rtol[i] = data->simulationInfo->tolerance;
    dasslData->atol[i] = data->simulationInfo->tolerance * fmax(fabs(data->modelData->realVarsData[i].attribute.nominal), 1e-32);
    infoStreamPrint(LOG_SOLVER, 0, "%d. %s -> %g", i+1, data->modelData->realVarsData[i].info.name, dasslData->atol[i]);
  }
  messageClose(LOG_SOLVER);


  /* let dassl return at every internal step */
  dasslData->info[2] = 1;


  /* define maximum step size, which is dassl is allowed to go */
  if (omc_flag[FLAG_MAX_STEP_SIZE])
  {
    double maxStepSize = atof(omc_flagValue[FLAG_MAX_STEP_SIZE]);

    assertStreamPrint(threadData, maxStepSize >= DASSL_STEP_EPS, "Selected maximum step size %e is too small.", maxStepSize);

    dasslData->rwork[1] = maxStepSize;
    dasslData->info[6] = 1;
    infoStreamPrint(LOG_SOLVER, 0, "maximum step size %g", dasslData->rwork[1]);
  }
  else
  {
    infoStreamPrint(LOG_SOLVER, 0, "maximum step size not set");
  }


  /* define initial step size, which is dassl is used every time it restarts */
  if (omc_flag[FLAG_INITIAL_STEP_SIZE])
  {
    double initialStepSize = atof(omc_flagValue[FLAG_INITIAL_STEP_SIZE]);

    assertStreamPrint(threadData, initialStepSize >= DASSL_STEP_EPS, "Selected initial step size %e is too small.", initialStepSize);

    dasslData->rwork[2] = initialStepSize;
    dasslData->info[7] = 1;
    infoStreamPrint(LOG_SOLVER, 0, "initial step size %g", dasslData->rwork[2]);
  }
  else
  {
    infoStreamPrint(LOG_SOLVER, 0, "initial step size not set");
  }


  /* define maximum integration order of dassl */
  if (omc_flag[FLAG_MAX_ORDER])
  {
    int maxOrder = atoi(omc_flagValue[FLAG_MAX_ORDER]);

    assertStreamPrint(threadData, maxOrder >= 1 && maxOrder <= 5, "Selected maximum order %d is out of range (1-5).", maxOrder);

    dasslData->iwork[2] = maxOrder;
    dasslData->info[8] = 1;
  }
  infoStreamPrint(LOG_SOLVER, 0, "maximum integration order %d", dasslData->info[8]?dasslData->iwork[2]:maxOrder);


  /* if FLAG_NOEQUIDISTANT_GRID is set, choose dassl step method */
  if (omc_flag[FLAG_NOEQUIDISTANT_GRID])
  {
    dasslData->dasslSteps = 1; /* TRUE */
    solverInfo->solverNoEquidistantGrid = 1;
  }
  else
  {
    dasslData->dasslSteps = 0; /* FALSE */
  }
  infoStreamPrint(LOG_SOLVER, 0, "use equidistant time grid %s", dasslData->dasslSteps?"NO":"YES");

  /* check if Flags FLAG_NOEQUIDISTANT_OUT_FREQ or FLAG_NOEQUIDISTANT_OUT_TIME are set */
  if (dasslData->dasslSteps){
    if (omc_flag[FLAG_NOEQUIDISTANT_OUT_FREQ])
    {
      dasslData->dasslStepsFreq = atoi(omc_flagValue[FLAG_NOEQUIDISTANT_OUT_FREQ]);
    }
    else if (omc_flag[FLAG_NOEQUIDISTANT_OUT_TIME])
    {
      dasslData->dasslStepsTime = atof(omc_flagValue[FLAG_NOEQUIDISTANT_OUT_TIME]);
      dasslData->rwork[1] = dasslData->dasslStepsTime;
      dasslData->info[6] = 1;
      infoStreamPrint(LOG_SOLVER, 0, "maximum step size %g", dasslData->rwork[1]);
    } else {
      dasslData->dasslStepsFreq = 1;
      dasslData->dasslStepsTime = 0.0;
    }

    if  (omc_flag[FLAG_NOEQUIDISTANT_OUT_FREQ] && omc_flag[FLAG_NOEQUIDISTANT_OUT_TIME]){
      warningStreamPrint(LOG_STDOUT, 0, "The flags are  \"noEquidistantOutputFrequency\" "
                                     "and \"noEquidistantOutputTime\" are in opposition "
                                     "to each other. The flag \"noEquidistantOutputFrequency\" superiors.");
     }
     infoStreamPrint(LOG_SOLVER, 0, "as the output frequency control is used: %d", dasslData->dasslStepsFreq);
     infoStreamPrint(LOG_SOLVER, 0, "as the output frequency time step control is used: %f", dasslData->dasslStepsTime);
  }

  /* if FLAG_DASSL_JACOBIAN is set, choose dassl jacobian calculation method */
  if (omc_flag[FLAG_DASSL_JACOBIAN])
  {
    for(i=1; i< DASSL_JAC_MAX;i++)
    {
      if(!strcmp((const char*)omc_flagValue[FLAG_DASSL_JACOBIAN], dasslJacobianMethodStr[i])){
        dasslData->dasslJacobian = (int)i;
        break;
      }
    }
    if(dasslData->dasslJacobian == DASSL_JAC_UNKNOWN)
    {
      if (ACTIVE_WARNING_STREAM(LOG_SOLVER))
      {
        warningStreamPrint(LOG_SOLVER, 1, "unrecognized jacobian calculation method %s, current options are:", (const char*)omc_flagValue[FLAG_DASSL_JACOBIAN]);
        for(i=1; i < DASSL_JAC_MAX; ++i)
        {
          warningStreamPrint(LOG_SOLVER, 0, "%-15s [%s]", dasslJacobianMethodStr[i], dasslJacobianMethodDescStr[i]);
        }
        messageClose(LOG_SOLVER);
      }
      throwStreamPrint(threadData,"unrecognized jacobian calculation method %s", (const char*)omc_flagValue[FLAG_DASSL_JACOBIAN]);
    }
  /* default case colored numerical jacobian */
  }
  else
  {
    dasslData->dasslJacobian = DASSL_COLOREDNUMJAC;
  }


  /* selects the calculation method of the jacobian */
  if(dasslData->dasslJacobian == DASSL_COLOREDNUMJAC ||
     dasslData->dasslJacobian == DASSL_COLOREDSYMJAC ||
     dasslData->dasslJacobian == DASSL_SYMJAC)
  {
    if (data->callback->initialAnalyticJacobianA(data, threadData))
    {
      infoStreamPrint(LOG_STDOUT, 0, "Jacobian or SparsePattern is not generated or failed to initialize! Switch back to normal.");
      dasslData->dasslJacobian = DASSL_INTERNALNUMJAC;
    }
  }
  /* default use a user sub-routine for JAC */
  dasslData->info[4] = 1;

  /* set up the appropriate function pointer */
  switch (dasslData->dasslJacobian){
    case DASSL_COLOREDNUMJAC:
      data->simulationInfo->jacobianEvals = data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern.maxColors;
      dasslData->jacobianFunction =  JacobianOwnNumColored;
      break;
    case DASSL_COLOREDSYMJAC:
      data->simulationInfo->jacobianEvals = data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern.maxColors;
      dasslData->jacobianFunction =  JacobianSymbolicColored;
      break;
    case DASSL_SYMJAC:
      dasslData->jacobianFunction =  JacobianSymbolic;
      break;
    case DASSL_NUMJAC:
      dasslData->jacobianFunction =  JacobianOwnNum;
      break;
    case DASSL_INTERNALNUMJAC:
      dasslData->jacobianFunction =  dummy_Jacobian;
      /* no user sub-routine for JAC */
      dasslData->info[4] = 0;
      break;
    default:
      throwStreamPrint(threadData,"unrecognized jacobian calculation method %s", (const char*)omc_flagValue[FLAG_DASSL_JACOBIAN]);
      break;
  }
  infoStreamPrint(LOG_SOLVER, 0, "jacobian is calculated by %s", dasslJacobianMethodDescStr[dasslData->dasslJacobian]);


  /* if FLAG_DASSL_NO_ROOTFINDING is set, choose dassl with out internal root finding */
  if(omc_flag[FLAG_DASSL_NO_ROOTFINDING])
  {
    dasslData->dasslRootFinding = 0;
    dasslData->zeroCrossingFunction = dummy_zeroCrossing;
    dasslData->ng = 0;
  }
  else
  {
    solverInfo->solverRootFinding = 1;
    dasslData->dasslRootFinding = 1;
    dasslData->zeroCrossingFunction = function_ZeroCrossingsDASSL;
  }
  infoStreamPrint(LOG_SOLVER, 0, "dassl uses internal root finding method %s", dasslData->dasslRootFinding?"YES":"NO");


  /* if FLAG_DASSL_NO_RESTART is set, choose dassl step method */
  if (omc_flag[FLAG_DASSL_NO_RESTART])
  {
    dasslData->dasslAvoidEventRestart = 1; /* TRUE */
  }
  else
  {
    dasslData->dasslAvoidEventRestart = 0; /* FALSE */
  }
  infoStreamPrint(LOG_SOLVER, 0, "dassl performs an restart after an event occurs %s", dasslData->dasslAvoidEventRestart?"NO":"YES");

  /* ### end configuration of dassl ### */


  messageClose(LOG_SOLVER);
  TRACE_POP
  return 0;
}


int dassl_deinitial(DASSL_DATA *dasslData)
{
  TRACE_PUSH
  unsigned int i;

  /* free work arrays for DASSL */
  free(dasslData->rwork);
  free(dasslData->iwork);
  free(dasslData->rpar);
  free(dasslData->ipar);
  free(dasslData->atol);
  free(dasslData->rtol);
  free(dasslData->info);
  free(dasslData->jroot);
  free(dasslData->ysave);
  free(dasslData->delta_hh);
  free(dasslData->newdelta);
  free(dasslData->stateDer);
  free(dasslData->dasslStatistics);
  free(dasslData->dasslStatisticsTmp);

  free(dasslData);

  TRACE_POP
  return 0;
}

/* \fn printCurrentStatesVector(int logLevel, double* y, DATA* data, double time)
 *
 * \param [in] [logLevel]
 * \param [in] [states]
 * \param [in] [data]
 * \param [in] [time]
 *
 * This function outputs states vector.
 *
 */
int printCurrentStatesVector(int logLevel, double* states, DATA* data, double time)
{
  int i;
  infoStreamPrint(logLevel, 1, "states at time=%g", time);
  for(i=0;i<data->modelData->nStates;++i)
  {
    infoStreamPrint(logLevel, 0, "%d. %s = %g", i+1, data->modelData->realVarsData[i].info.name, states[i]);
  }
  messageClose(logLevel);

  return 0;
}



/**********************************************************************************************
 * DASSL with synchronous treating of when equation
 *   - without integrated ZeroCrossing method.
 *   + ZeroCrossing are handled outside DASSL.
 *   + if no event occurs outside DASSL performs a warm-start
 **********************************************************************************************/
int dassl_step(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo)
{
  TRACE_PUSH
  double tout = 0;
  int i = 0;
  unsigned int ui = 0;
  int retVal = 0;
  int saveJumpState;
  static unsigned int dasslStepsOutputCounter = 1;

  DASSL_DATA *dasslData = (DASSL_DATA*) solverInfo->solverData;

  SIMULATION_DATA *sData = data->localData[0];
  SIMULATION_DATA *sDataOld = data->localData[1];

  MODEL_DATA *mData = (MODEL_DATA*) data->modelData;

  modelica_real* stateDer = dasslData->stateDer;

  dasslData->rpar[0] = (double*) (void*) data;
  dasslData->rpar[1] = (double*) (void*) dasslData;
  dasslData->rpar[2] = (double*) (void*) threadData;


//  printf("\tlastdesiredStep: %f\n", solverInfo->lastdesiredStep);
//  printf("\tstateEvents: %d\n", solverInfo->stateEvents);
//  printf("\tdidEventStep: %d\n", solverInfo->didEventStep);
//  printf("\tsampleEvents: %d\n", solverInfo->sampleEvents);
//  printf("\tintegratorSteps: %d\n", solverInfo->integratorSteps);

  saveJumpState = threadData->currentErrorStage;
  threadData->currentErrorStage = ERROR_INTEGRATOR;

  /* try */
#if !defined(OMC_EMCC)
  MMC_TRY_INTERNAL(simulationJumpBuffer)
#endif

  assertStreamPrint(threadData, 0 != dasslData->rpar, "could not passed to DDASKR");

  /* If an event is triggered and processed restart dassl. */
  if(!dasslData->dasslAvoidEventRestart && (solverInfo->didEventStep || 0 == dasslData->idid))
  {
    debugStreamPrint(LOG_EVENTS_V, 0, "Event-management forced reset of DDASKR");
    /* obtain reset */
    dasslData->info[0] = 0;
    dasslData->idid = 0;

    memcpy(stateDer, data->localData[1]->realVars + data->modelData->nStates, sizeof(double)*data->modelData->nStates);
  }

  /* Calculate steps until TOUT is reached */
  if (dasslData->dasslSteps)
  {
    /* If dasslsteps is selected, the dassl run to stopTime or next sample event */
    if (data->simulationInfo->nextSampleEvent < data->simulationInfo->stopTime)
    {
      tout = data->simulationInfo->nextSampleEvent;
    }
    else
    {
      tout = data->simulationInfo->stopTime;
    }
  }
  else
  {
    tout = solverInfo->currentTime + solverInfo->currentStepSize;
  }

  /* Check that tout is not less than timeValue
   * else will dassl get in trouble. If that is the case we skip the current step. */
  if (solverInfo->currentStepSize < DASSL_STEP_EPS)
  {
    infoStreamPrint(LOG_DASSL, 0, "Desired step to small try next one");
    infoStreamPrint(LOG_DASSL, 0, "Interpolate linear");

    /*euler step*/
    for(i = 0; i < data->modelData->nStates; i++)
    {
      sData->realVars[i] = sDataOld->realVars[i] + stateDer[i] * solverInfo->currentStepSize;
    }
    sData->timeValue = solverInfo->currentTime + solverInfo->currentStepSize;
    data->callback->functionODE(data, threadData);
    solverInfo->currentTime = sData->timeValue;

    TRACE_POP
    return 0;
  }

  do
  {
    infoStreamPrint(LOG_DASSL, 1, "new step at time = %.15g", solverInfo->currentTime);

    /* rhs final flag is FALSE during for dassl evaluation */
    RHSFinalFlag = 0;

    /* read input vars */
    externalInputUpdate(data);
    data->callback->input_function(data, threadData);

    DDASKR(functionODE_residual, (int*) &mData->nStates,
            &solverInfo->currentTime, sData->realVars, stateDer, &tout,
            dasslData->info, dasslData->rtol, dasslData->atol, &dasslData->idid,
            dasslData->rwork, &dasslData->lrw, dasslData->iwork, &dasslData->liw,
            (double*) (void*) dasslData->rpar, dasslData->ipar, dasslData->jacobianFunction, dummy_precondition,
            dasslData->zeroCrossingFunction, (int*) &dasslData->ng, dasslData->jroot);

    /* closing new step message */
    messageClose(LOG_DASSL);

    /* set ringbuffer time to current time */
    sData->timeValue = solverInfo->currentTime;

    /* rhs final flag is TRUE during for output evaluation */
    RHSFinalFlag = 1;

    if(dasslData->idid == -1)
    {
      fflush(stderr);
      fflush(stdout);
      warningStreamPrint(LOG_DASSL, 0, "A large amount of work has been expended.(About 500 steps). Trying to continue ...");
      infoStreamPrint(LOG_DASSL, 0, "DASSL will try again...");
      dasslData->info[0] = 1; /* try again */
      if (solverInfo->currentTime <= data->simulationInfo->stopTime)
        continue;
    }
    else if(dasslData->idid < 0)
    {
      fflush(stderr);
      fflush(stdout);
      retVal = continue_DASSL(&dasslData->idid, &data->simulationInfo->tolerance);
      warningStreamPrint(LOG_STDOUT, 0, "can't continue. time = %f", sData->timeValue);
      TRACE_POP
      break;
    }
    else if(dasslData->idid == 5)
    {
      threadData->currentErrorStage = ERROR_EVENTSEARCH;
    }

    /* emit step, if dasslsteps is selected */
    if (dasslData->dasslSteps)
    {
      if (omc_flag[FLAG_NOEQUIDISTANT_OUT_FREQ]){
        /* output every n-th time step */
        if (dasslStepsOutputCounter >= dasslData->dasslStepsFreq){
          dasslStepsOutputCounter = 1; /* next line set it to one */
          break;
        }
        dasslStepsOutputCounter++;
      } else if (omc_flag[FLAG_NOEQUIDISTANT_OUT_TIME]){
        /* output when time>=k*timeValue */
        if (solverInfo->currentTime > dasslStepsOutputCounter * dasslData->dasslStepsTime){
          dasslStepsOutputCounter++;
          break;
        }
      } else {
        break;
      }
    }

  } while(dasslData->idid == 1);

#if !defined(OMC_EMCC)
  MMC_CATCH_INTERNAL(simulationJumpBuffer)
#endif
  threadData->currentErrorStage = saveJumpState;

  /* if a state event occurs than no sample event does need to be activated  */
  if (data->simulationInfo->sampleActivated && solverInfo->currentTime < data->simulationInfo->nextSampleEvent)
  {
    data->simulationInfo->sampleActivated = 0;
  }


  if(ACTIVE_STREAM(LOG_DASSL))
  {
    infoStreamPrint(LOG_DASSL, 1, "dassl call statistics: ");
    infoStreamPrint(LOG_DASSL, 0, "value of idid: %d", (int)dasslData->idid);
    infoStreamPrint(LOG_DASSL, 0, "current time value: %0.4g", solverInfo->currentTime);
    infoStreamPrint(LOG_DASSL, 0, "current integration time value: %0.4g", dasslData->rwork[3]);
    infoStreamPrint(LOG_DASSL, 0, "step size H to be attempted on next step: %0.4g", dasslData->rwork[2]);
    infoStreamPrint(LOG_DASSL, 0, "step size used on last successful step: %0.4g", dasslData->rwork[6]);
    infoStreamPrint(LOG_DASSL, 0, "the order of the method used on the last step: %d", dasslData->iwork[7]);
    infoStreamPrint(LOG_DASSL, 0, "the order of the method to be attempted on the next step: %d", dasslData->iwork[8]);
    infoStreamPrint(LOG_DASSL, 0, "number of steps taken so far: %d", (int)dasslData->iwork[10]);
    infoStreamPrint(LOG_DASSL, 0, "number of calls of functionODE() : %d", (int)dasslData->iwork[11]);
    infoStreamPrint(LOG_DASSL, 0, "number of calculation of jacobian : %d", (int)dasslData->iwork[12]);
    infoStreamPrint(LOG_DASSL, 0, "total number of convergence test failures: %d", (int)dasslData->iwork[13]);
    infoStreamPrint(LOG_DASSL, 0, "total number of error test failures: %d", (int)dasslData->iwork[14]);
    messageClose(LOG_DASSL);
  }

  /* save dassl stats */
  for(ui = 0; ui < numStatistics; ui++)
  {
    assert(10 + ui < dasslData->liw);
    dasslData->dasslStatisticsTmp[ui] = dasslData->iwork[10 + ui];
  }

  infoStreamPrint(LOG_DASSL, 0, "Finished DASSL step.");

  TRACE_POP
  return retVal;
}

static int
continue_DASSL(int* idid, double* atol)
{
  TRACE_PUSH
  int retValue = -1;

  switch(*idid)
  {
  case 1:
  case 2:
  case 3:
    /* 1-4 means success */
    break;
  case -1:
    warningStreamPrint(LOG_DASSL, 0, "A large amount of work has been expended.(About 500 steps). Trying to continue ...");
    retValue = 1; /* adrpo: try to continue */
    break;
  case -2:
    warningStreamPrint(LOG_STDOUT, 0, "The error tolerances are too stringent");
    retValue = -2;
    break;
  case -3:
    /* wbraun: don't throw at this point let the solver handle it */
    /* throwStreamPrint("DDASKR: THE LAST STEP TERMINATED WITH A NEGATIVE IDID value"); */
    retValue = -3;
    break;
  case -6:
    warningStreamPrint(LOG_STDOUT, 0, "DDASSL had repeated error test failures on the last attempted step.");
    retValue = -6;
    break;
  case -7:
    warningStreamPrint(LOG_STDOUT, 0, "The corrector could not converge.");
    retValue = -7;
    break;
  case -8:
    warningStreamPrint(LOG_STDOUT, 0, "The matrix of partial derivatives is singular.");
    retValue = -8;
    break;
  case -9:
    warningStreamPrint(LOG_STDOUT, 0, "The corrector could not converge. There were repeated error test failures in this step.");
    retValue = -9;
    break;
  case -10:
    warningStreamPrint(LOG_STDOUT, 0, "A Modelica assert prevents the integrator to continue. For more information use -lv LOG_SOLVER");
    retValue = -10;
    break;
  case -11:
    warningStreamPrint(LOG_STDOUT, 0, "IRES equal to -2 was encountered and control is being returned to the calling program.");
    retValue = -11;
    break;
  case -12:
    warningStreamPrint(LOG_STDOUT, 0, "DDASSL failed to compute the initial YPRIME.");
    retValue = -12;
    break;
  case -33:
    warningStreamPrint(LOG_STDOUT, 0, "The code has encountered trouble from which it cannot recover.");
    retValue = -33;
    break;
  }

  TRACE_POP
  return retValue;
}


int functionODE_residual(double *t, double *y, double *yd, double* cj, double *delta,
                    int *ires, double *rpar, int *ipar)
{
  TRACE_PUSH
  DATA* data = (DATA*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)((double**)rpar)[2];

  double timeBackup;
  long i;
  int saveJumpState;
  int success = 0;

  if (data->simulationInfo->currentContext == CONTEXT_ALGEBRAIC)
  {
    setContext(data, t, CONTEXT_ODE);
  }
  printCurrentStatesVector(LOG_DASSL_STATES, y, data, *t);

  timeBackup = data->localData[0]->timeValue;
  data->localData[0]->timeValue = *t;

  saveJumpState = threadData->currentErrorStage;
  threadData->currentErrorStage = ERROR_INTEGRATOR;

  /* try */
#if !defined(OMC_EMCC)
  MMC_TRY_INTERNAL(simulationJumpBuffer)
#endif

  /* read input vars */
  externalInputUpdate(data);
  data->callback->input_function(data, threadData);

  /* eval input vars */
  data->callback->functionODE(data, threadData);

  /* get the difference between the temp_xd(=localData->statesDerivatives)
     and xd(=statesDerivativesBackup) */
  for(i=0; i < data->modelData->nStates; i++)
  {
    delta[i] = data->localData[0]->realVars[data->modelData->nStates + i] - yd[i];
  }
  success = 1;
#if !defined(OMC_EMCC)
  MMC_CATCH_INTERNAL(simulationJumpBuffer)
#endif

  if (!success) {
    *ires = -1;
  }

  threadData->currentErrorStage = saveJumpState;

  data->localData[0]->timeValue = timeBackup;

  if (data->simulationInfo->currentContext == CONTEXT_ODE){
    unsetContext(data);
  }

  TRACE_POP
  return 0;
}

int function_ZeroCrossingsDASSL(int *neqm, double *t, double *y, double *yp,
        int *ng, double *gout, double *rpar, int* ipar)
{
  TRACE_PUSH
  DATA* data = (DATA*)(void*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)(void*)((double**)rpar)[2];

  double timeBackup;
  int saveJumpState;

  if (data->simulationInfo->currentContext == CONTEXT_ALGEBRAIC)
  {
    setContext(data, t, CONTEXT_EVENTS);
  }

  saveJumpState = threadData->currentErrorStage;
  threadData->currentErrorStage = ERROR_EVENTSEARCH;

  timeBackup = data->localData[0]->timeValue;
  data->localData[0]->timeValue = *t;

  /* read input vars */
  externalInputUpdate(data);
  data->callback->input_function(data, threadData);
  /* eval needed equations*/
  data->callback->function_ZeroCrossingsEquations(data, threadData);

  data->callback->function_ZeroCrossings(data, threadData, gout);

  threadData->currentErrorStage = saveJumpState;
  data->localData[0]->timeValue = timeBackup;

  if (data->simulationInfo->currentContext == CONTEXT_EVENTS){
    unsetContext(data);
  }

  TRACE_POP
  return 0;
}


int functionJacAColored(DATA* data, threadData_t *threadData, double* jac)
{
  TRACE_PUSH
  const int index = data->callback->INDEX_JAC_A;
  unsigned int i,j,l,k,ii;

  for(i=0; i < data->simulationInfo->analyticJacobians[index].sparsePattern.maxColors; i++)
  {
    for(ii=0; ii < data->simulationInfo->analyticJacobians[index].sizeCols; ii++)
      if(data->simulationInfo->analyticJacobians[index].sparsePattern.colorCols[ii]-1 == i)
        data->simulationInfo->analyticJacobians[index].seedVars[ii] = 1;

    data->callback->functionJacA_column(data, threadData);

    for(j = 0; j < data->simulationInfo->analyticJacobians[index].sizeCols; j++)
    {
      if(data->simulationInfo->analyticJacobians[index].seedVars[j] == 1)
      {
        if(j==0)
          ii = 0;
        else
          ii = data->simulationInfo->analyticJacobians[index].sparsePattern.leadindex[j-1];
        while(ii < data->simulationInfo->analyticJacobians[index].sparsePattern.leadindex[j])
        {
          l  = data->simulationInfo->analyticJacobians[index].sparsePattern.index[ii];
          k  = j*data->simulationInfo->analyticJacobians[index].sizeRows + l;
          jac[k] = data->simulationInfo->analyticJacobians[index].resultVars[l];
          ii++;
        };
      }
    }
    for(ii=0; ii < data->simulationInfo->analyticJacobians[index].sizeCols; ii++)
      if(data->simulationInfo->analyticJacobians[index].sparsePattern.colorCols[ii]-1 == i) data->simulationInfo->analyticJacobians[index].seedVars[ii] = 0;

  }

  TRACE_POP
  return 0;
}


int functionJacASym(DATA* data, threadData_t *threadData, double* jac)
{
  TRACE_PUSH
  const int index = data->callback->INDEX_JAC_A;
  unsigned int i,j,k;

  k = 0;
  for(i=0; i < data->simulationInfo->analyticJacobians[index].sizeCols; i++)
  {
    data->simulationInfo->analyticJacobians[index].seedVars[i] = 1.0;

    data->callback->functionJacA_column(data, threadData);

    for(j = 0; j < data->simulationInfo->analyticJacobians[index].sizeRows; j++)
    {
      jac[k++] = data->simulationInfo->analyticJacobians[index].resultVars[j];
    }

    data->simulationInfo->analyticJacobians[index].seedVars[i] = 0.0;
  }

  TRACE_POP
  return 0;
}

/*
 * provides a analytical Jacobian to be used with DASSL
 */

static int JacobianSymbolicColored(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
         double *rpar, int* ipar)
{
  TRACE_PUSH
  DATA* data = (DATA*)(void*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)(void*)((double**)rpar)[2];
  double* backupStates;
  double timeBackup;
  int i;
  int j;

  setContext(data, t, CONTEXT_JACOBIAN);

  backupStates = data->localData[0]->realVars;
  timeBackup = data->localData[0]->timeValue;

  data->localData[0]->timeValue = *t;
  data->localData[0]->realVars = y;
  /* read input vars */
  externalInputUpdate(data);
  data->callback->input_function(data, threadData);
  /* eval ode*/
  data->callback->functionODE(data, threadData);
  functionJacAColored(data, threadData, pd);

  /* add cj to the diagonal elements of the matrix */
  j = 0;
  for(i = 0; i < data->modelData->nStates; i++)
  {
    pd[j] -= (double) *cj;
    j += data->modelData->nStates + 1;
  }
  data->localData[0]->realVars = backupStates;
  data->localData[0]->timeValue = timeBackup;

  unsetContext(data);

  TRACE_POP
  return 0;
}

/*
 * provides a analytical Jacobian to be used with DASSL
 */
static int JacobianSymbolic(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
         double *rpar, int* ipar)
{
  TRACE_PUSH
  DATA* data = (DATA*)(void*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)(void*)((double**)rpar)[2];
  double* backupStates;
  double timeBackup;
  int i;
  int j;

  setContext(data, t, CONTEXT_JACOBIAN);

  backupStates = data->localData[0]->realVars;
  timeBackup = data->localData[0]->timeValue;

  data->localData[0]->timeValue = *t;
  data->localData[0]->realVars = y;
  /* read input vars */
  externalInputUpdate(data);
  data->callback->input_function(data, threadData);
  /* eval ode*/
  data->callback->functionODE(data, threadData);
  functionJacASym(data, threadData, pd);

  /* add cj to the diagonal elements of the matrix */
  j = 0;
  for(i = 0; i < data->modelData->nStates; i++)
  {
    pd[j] -= (double) *cj;
    j += data->modelData->nStates + 1;
  }
  data->localData[0]->realVars = backupStates;
  data->localData[0]->timeValue = timeBackup;

  unsetContext(data);

  TRACE_POP
  return 0;
}

/*
 *  function calculates a jacobian matrix by
 *  numerical method finite differences
 */
int jacA_num(DATA* data, double *t, double *y, double *yprime, double *delta, double *matrixA, double *cj, double *h, double *wt, double *rpar, int *ipar)
{
  TRACE_PUSH
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];

  double delta_h = dasslData->sqrteps;
  double delta_hh,delta_hhh, deltaInv;
  double ysave;
  int ires;
  int i,j;


  for(i=data->modelData->nStates-1; i >= 0; i--)
  {
    delta_hhh = *h * yprime[i];
    delta_hh = delta_h * fmax(fmax(fabs(y[i]),fabs(delta_hhh)),fabs(1. / wt[i]));
    delta_hh = (delta_hhh >= 0 ? delta_hh : -delta_hh);
    delta_hh = y[i] + delta_hh - y[i];
    deltaInv = 1. / delta_hh;
    ysave = y[i];
    y[i] += delta_hh;

    /* internal dassl numerical jacobian is
     * calculated by adding cj to yprime.
     * This lead to numerical cancellations.
     */
    /*yprime[i] += *cj * delta_hh;*/

    functionODE_residual(t, y, yprime, cj, dasslData->newdelta, &ires, rpar, ipar);

    for(j = data->modelData->nStates-1; j >= 0 ; j--)
    {
      matrixA[i*data->modelData->nStates+j] = (dasslData->newdelta[j] - delta[j]) * deltaInv;
    }
    y[i] = ysave;
  }

  TRACE_POP
  return 0;
}

/*
 * provides a numerical Jacobian to be used with DASSL
 */
static int JacobianOwnNum(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
    double *rpar, int* ipar)
{
  TRACE_PUSH
  int i,j;
  DATA* data = (DATA*)(void*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)(void*)((double**)rpar)[2];

  setContext(data, t, CONTEXT_JACOBIAN);

  if(jacA_num(data, t, y, yprime, deltaD, pd, cj, h, wt, rpar, ipar))
  {
    throwStreamPrint(threadData, "Error, can not get Matrix A ");
    TRACE_POP
    return 1;
  }
  j = 0;
  /* add cj to diagonal elements and store in pd */
  for(i = 0; i < data->modelData->nStates; i++)
  {
    pd[j] -= (double) *cj;
    j += data->modelData->nStates + 1;
  }
  unsetContext(data);

  TRACE_POP
  return 0;
}


/*
 *  function calculates a jacobian matrix by
 *  numerical method finite differences
 */
int jacA_numColored(DATA* data, double *t, double *y, double *yprime, double *delta, double *matrixA, double *cj, double *h, double *wt, double *rpar, int *ipar)
{
  TRACE_PUSH
  const int index = data->callback->INDEX_JAC_A;
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  double delta_h = dasslData->sqrteps;
  double delta_hhh;
  int ires;
  double* delta_hh = dasslData->delta_hh;
  double* ysave = dasslData->ysave;

  unsigned int i,j,l,k,ii;

  for(i = 0; i < data->simulationInfo->analyticJacobians[index].sparsePattern.maxColors; i++)
  {
    for(ii=0; ii < data->simulationInfo->analyticJacobians[index].sizeCols; ii++)
    {
      if(data->simulationInfo->analyticJacobians[index].sparsePattern.colorCols[ii]-1 == i)
      {
        delta_hhh = *h * yprime[ii];
        delta_hh[ii] = delta_h * fmax(fmax(fabs(y[ii]),fabs(delta_hhh)),fabs(1./wt[ii]));
        delta_hh[ii] = (delta_hhh >= 0 ? delta_hh[ii] : -delta_hh[ii]);
        delta_hh[ii] = y[ii] + delta_hh[ii] - y[ii];

        ysave[ii] = y[ii];
        y[ii] += delta_hh[ii];

        delta_hh[ii] = 1. / delta_hh[ii];
      }
    }

    functionODE_residual(t, y, yprime, cj, dasslData->newdelta, &ires, rpar, ipar);

    increaseJacContext(data);

    for(ii = 0; ii < data->simulationInfo->analyticJacobians[index].sizeCols; ii++)
    {
      if(data->simulationInfo->analyticJacobians[index].sparsePattern.colorCols[ii]-1 == i)
      {
        if(ii==0)
          j = 0;
        else
          j = data->simulationInfo->analyticJacobians[index].sparsePattern.leadindex[ii-1];
        while(j < data->simulationInfo->analyticJacobians[index].sparsePattern.leadindex[ii])
        {
          l  =  data->simulationInfo->analyticJacobians[index].sparsePattern.index[j];
          k  = l + ii*data->simulationInfo->analyticJacobians[index].sizeRows;
          matrixA[k] = (dasslData->newdelta[l] - delta[l]) * delta_hh[ii];
          j++;
        };
        y[ii] = ysave[ii];
      }
    }
  }

  TRACE_POP
  return 0;
}

/*
 * provides a numerical Jacobian to be used with DASSL
 */
static int JacobianOwnNumColored(double *t, double *y, double *yprime, double *deltaD, double *pd, double *cj, double *h, double *wt,
   double *rpar, int* ipar)
{
  TRACE_PUSH
  DATA* data = (DATA*)(void*)((double**)rpar)[0];
  DASSL_DATA* dasslData = (DASSL_DATA*)(void*)((double**)rpar)[1];
  threadData_t *threadData = (threadData_t*)(void*)((double**)rpar)[2];
  int i,j;

  setContext(data, t, CONTEXT_JACOBIAN);

  if(jacA_numColored(data, t, y, yprime, deltaD, pd, cj, h, wt, rpar, ipar))
  {
    throwStreamPrint(threadData, "Error, can not get Matrix A ");
    TRACE_POP
    return 1;
  }

  /* add cj to diagonal elements and store in pd */
  j = 0;
  for(i = 0; i < data->modelData->nStates; i++)
  {
    pd[j] -= (double) *cj;
    j += data->modelData->nStates + 1;
  }
  unsetContext(data);

  TRACE_POP
  return 0;
}

#ifdef __cplusplus
}
#endif
