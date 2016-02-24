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

#include <stdio.h>
#include "solver_main.h"

#include "simulation/simulation_runtime.h"
#include "simulation/results/simulation_result.h"
#include "openmodelica_func.h"

#include "util/omc_error.h"
#include "simulation/options.h"


/*! enum error_msg
 * \brief  Returnvalues of the functions in this file
 */
enum error_msg
{
  ISNAN = -3L,      /*!< Time of next change is #QNAN. */
  UNKNOWN = -2L,    /*!< Unspecific error. */
  OO_MEMORY = -1L,  /*!< Allocation of memory fails. */
  OK = 0L           /*!< Everything is fine. */
};

const modelica_real EPS = 1e-15;
const modelica_real deltaQFactor =0.0001;


/* Needed if we want to write all the variables into a file*/
/* #define D */

static modelica_integer deltaQ( DATA* data,const modelica_real dQ, const modelica_integer index, modelica_real* dTnextQ, modelica_real* nextQ, modelica_real* diffQ);
static modelica_integer getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k);
static modelica_integer getStatesInDer(const unsigned int* index, const unsigned int* leadindex, const uinteger ROWS, const uinteger STATES, uinteger** StatesInDer);
static modelica_integer qss_step(DATA* data, SOLVER_INFO* solverInfo);
static uinteger minStep( const modelica_real* tqp, const uinteger size );

/*! performQSSSimulation(DATA* data, SOLVER_INFO* solverInfo)
 *
 *  \param [ref] [data]
 *  \param [ref] [solverInfo]
 *
 *  This function performs the simulation controlled by solverInfo.
 */
modelica_integer prefixedName_performQSSSimulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo)
{
  TRACE_PUSH
  printf("liqss.c: First comment in performQSSSimulation\n\n");
  SIMULATION_INFO *simInfo = data->simulationInfo;

  MODEL_DATA *mData = data->modelData;
  uinteger currStepNo = 0;
  modelica_integer retValIntegrator = 0;
  modelica_integer retValue = 0;
  uinteger ind = 0;

  solverInfo->currentTime = simInfo->startTime;
  warningStreamPrint(LOG_STDOUT, 0, "This QSS method is under development and should not be used yet.");

  if (data->callback->initialAnalyticJacobianA(data, threadData))
  {
    infoStreamPrint(LOG_STDOUT, 0, "Jacobian or sparse pattern is not generated or failed to initialize.");
    return UNKNOWN;
  }
  //printSparseStructure(data, LOG_SOLVER);

/* *********************************************************************************** */
  /* Initialization */
  uinteger i = 0; /* loop var */
  SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
  modelica_real* state = sData->realVars;
  modelica_real* stateDer = sData->realVars + data->modelData->nStates;
  // Sit my eie shit hier in
  //modelica_real* antonDer = data->modelData->nStates;



  const SPARSE_PATTERN* pattern = &(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern);
  const uinteger ROWS = data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sizeRows;
  const uinteger STATES = data->modelData->nStates;
  uinteger numDer = 0;  /* number of derivatives influenced by state k */

  modelica_boolean fail = 0;
  modelica_real* qik = NULL;  /* Approximation of states */
  modelica_real* xik = NULL;  /* states */
  modelica_real* derXik = NULL;  /* Derivative of states */
  modelica_real* tq = NULL;    /* Time of approximation, because not all approximations are calculated at a specific time, each approx. has its own timestamp */
  modelica_real* tx = NULL;    /* Time of the states, because not all states are calculated at a specific time, each state has its own timestamp */
  modelica_real* tqp = NULL;    /* Time of the next change in state */
  modelica_real* nQh = NULL;    /* next value of the state */
  modelica_real* dQ = NULL;    /* change in quantity of every state, default = nominal*10^-4 */
//   printf("---------------###---------------\n\n");
//   printf("Some stats:\n\n");
//   printf("(A) Firstly, we consider some info of the object DATA *data:\n");
//   printf("(1) SIMULATION_INFO (*simulationInfo):\n");
//   printf("\tstartTime: %f\n", data->simulationInfo->startTime);
//   printf("\tstopTime: %f\n", data->simulationInfo->stopTime);
//   printf("\tnumSteps: %d\n", data->simulationInfo->numSteps);
//   printf("\ttolerance: %f\n", data->simulationInfo->tolerance);
//   printf("\ttolerance: %f\n", (*(*data).simulationInfo).tolerance);
//   printf("\t*solverMethod: %c\n", *(data->simulationInfo->solverMethod));
//   printf("\t*outputFormat: %c\n", *(data->simulationInfo->outputFormat));
//   printf("\t*variableFilter: %c\n", *(data->simulationInfo->variableFilter));
//   printf("\tlsMethod: %d\n", data->simulationInfo->lsMethod);
//   printf("\tmixedMethod: %d\n", data->simulationInfo->mixedMethod);
//   printf("\tnlsMethod: %d\n", data->simulationInfo->nlsMethod);
//   printf("\tnewtonStrategy: %d\n", data->simulationInfo->newtonStrategy);
//   printf("\tnlsCsvInfomation: %d\n", data->simulationInfo->nlsCsvInfomation);
//   printf("\tcurrentContext: %d\n", data->simulationInfo->currentContext);
//   printf("\tcurrentContextOld: %d\n", data->simulationInfo->currentContextOld);
//   printf("\tjacobianEvals: %d\n", data->simulationInfo->jacobianEvals);
//   printf("\tlambda: %d\n", data->simulationInfo->lambda);
//   printf("\tnextSampleEvent: %f\n", data->simulationInfo->nextSampleEvent);
//   printf("\tnextSampleTimes[0]: %f\n", data->simulationInfo->nextSampleTimes[0]);
//   printf("\tnextSampleTimes[1]: %f\n", data->simulationInfo->nextSampleTimes[1]);
//   printf("\tnextSampleTimes[2]: %f\n", data->simulationInfo->nextSampleTimes[2]);
//
//
//   printf("\n\tThere are even more although I believe they will not be useful for now.\n\n");
//
//   printf("(2) SIMULATION_DATA (**localData):\n");
//   printf("\ttimeValue: %f\n", data->localData[0]->timeValue);
//   printf("The following are the initial values.\n");
//   printf("\trealVars[0]: %f\n", data->localData[0]->realVars[0]);
//   printf("\trealVars[1]: %f\n", data->localData[0]->realVars[1]);
//   printf("\tintegerVars: %d\n", data->localData[0]->integerVars[1]);
//   printf("\t*booleanVars: %s\n", *(data->localData[0]->booleanVars));
//   printf("\t*stringVars: %s\n", *(data->localData[0]->stringVars));
//   //printf("\t*inlineVars: %f\n", *(data->localData[0]->inlineVars)); //Errors when attempting to print this one. Not sure why though.
//   printf("\n\tI am not sure this garbage struct is necessary. One \'timeValue\' is used. I am sure we can use a local variable or a similar route.\n\n");
//
//   printf("(3) MODEL_DATA (*modelData):\n");
//   printf("\trealVarsData[0].attribute.nominal: %f\n", data->modelData->realVarsData[0].attribute.nominal);
//   printf("\trealVarsData[1].attribute.nominal: %f\n", data->modelData->realVarsData[1].attribute.nominal);
//   printf("\trealVarsData[0].attribute.useStart: %s\n", data->modelData->realVarsData[0].attribute.useStart ? "true" : "false");
//   printf("\trealVarsData[1].attribute.useStart: %s\n", data->modelData->realVarsData[1].attribute.useStart ? "true" : "false");
//   printf("\tnStates: %d\n", data->modelData->nStates);
//   printf("\tnVariablesReal: %d\n", data->modelData->nVariablesReal);
//   printf("\tnDiscreteReal: %d\n", data->modelData->nDiscreteReal);
//   printf("\tnVariablesInteger: %d\n", data->modelData->nVariablesInteger);
//   printf("\tnVariablesBoolean: %d\n", data->modelData->nVariablesBoolean);
//   printf("\tnVariablesString: %d\n", data->modelData->nVariablesString);
//   printf("\tnParametersReal: %d\n", data->modelData->nParametersReal);
//   printf("\tnVariablesReal: %d\n", data->modelData->nVariablesReal);
//   printf("\tnParametersInteger: %d\n", data->modelData->nParametersInteger);
//   printf("\tnParametersBoolean: %d\n", data->modelData->nParametersBoolean);
//   printf("\tnParametersString: %d\n", data->modelData->nParametersString);
//   printf("\tnInputVars: %d\n", data->modelData->nInputVars);
//   printf("\tnOutputVars: %d\n", data->modelData->nOutputVars);
//   printf("\tnZeroCrossings: %d\n", data->modelData->nZeroCrossings);
//   printf("\tnMathEvents: %d\n", data->modelData->nMathEvents);
//   printf("\tnDelayExpressions: %d\n", data->modelData->nDelayExpressions);
//   printf("\tnExtObjs: %d\n", data->modelData->nExtObjs);
//   printf("\tnMixedSystems: %d\n", data->modelData->nMixedSystems);
//   printf("\tnLinearSystems: %d\n", data->modelData->nLinearSystems);
//   printf("\tnNonLinearSystems: %d\n", data->modelData->nNonLinearSystems);
//   printf("\tnStateSets: %d\n", data->modelData->nStateSets);
//   printf("\tnInlineVars: %d\n", data->modelData->nInlineVars);
//   printf("\tnOptimizeConstraints: %d\n", data->modelData->nOptimizeConstraints);
//   printf("\tnOptimizeFinalConstraints: %d\n\n", data->modelData->nOptimizeFinalConstraints);
//
//   printf("(4) RINGBUFFER* (simulationData) -- The first comment in the .c file is that this does not work. This struct is never used in the qss solver.\n\n");
//   //printf("\titemSize: %d\n", data->simulationData->itemSize);
//
//   printf("(5) OpenModelicaGeneratedFunctionCallbacks (*callback). Just a number of functions updating the results. Not even sure where all these functions are defined. We may want to investigate where these functions are defined if we want to have a full understanding of the how OpenModelica enables it's solvers.\n\n");
//
//   printf("(B) threadData_t threadData. Next we have to consider the object threadData_t. However, I cannot find the definition of this object. Must investigate and get all attributes.\n\n");
//
//   printf("(C) SOLVER_INFO solverInfo:\n");
//   printf("\tcurrentTime: %f\n", solverInfo->currentTime);
//   printf("\tcurrentStepSize: %f\n", solverInfo->currentStepSize);
//   printf("\tlaststep: %f\n", solverInfo->laststep);
//   printf("\tsolverMethod: %d\n", solverInfo->solverMethod);
//   printf("\tsolverStepSize: %f\n", solverInfo->solverStepSize);
//   printf("\tsolverRootFinding: %s\n", solverInfo->solverRootFinding ? "true" : "false");
//   printf("\tsolverNoEquidistantGrid: %s\n", solverInfo->solverNoEquidistantGrid ? "true" : "false");
//   printf("\tlastdesiredStep: %f\n", solverInfo->lastdesiredStep);
//   printf("\tstateEvents: %d\n", solverInfo->stateEvents);
//   printf("\tdidEventStep: %d\n", solverInfo->didEventStep);
// //  printf("\teventLst.itemSize %d\n", solverInfo->eventLst->itemSize);
// //  printf("\teventLst.length %d\n", solverInfo->eventLst.length);
//   printf("\tsampleEvents: %d\n", solverInfo->sampleEvents);
//   printf("\tintegratorSteps: %d\n", solverInfo->integratorSteps);
//
//
   printf("---------------###---------------\n\n");
   for (i = 0; i < STATES; i++){
	   printf("state[i]: %d %f\n",i,state[i]);
   }
   for (i = 0; i < STATES; i++){
   	   printf("stateDer[i]: %d %f\n",i,stateDer[i]);
      }
   for (i = 0; i < STATES; i++){
	   state[i]= state[i] +stateDer[i]*0.45;
   }

   solverInfo->currentTime = 0.45;
   /* update continous system */
       sData->timeValue = solverInfo->currentTime;
       externalInputUpdate(data);
       data->callback->input_function(data, threadData);
       data->callback->functionODE(data, threadData);
       data->callback->functionAlgebraics(data, threadData);
       data->callback->output_function(data, threadData);
       data->callback->function_storeDelayed(data, threadData);
       data->callback->functionDAE(data, threadData);
   printf("---------------###---------------\n\n");
      for (i = 0; i < STATES; i++){
   	   printf("state[i]: %d %f\n",i,state[i]);
      }
      for (i = 0; i < STATES; i++){
      	   printf("stateDer[i]: %d %f\n",i,stateDer[i]);
         }


  /* allocate memory*/
  qik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (qik == NULL) ? 1 : ( 0 | fail);
  xik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (xik == NULL) ? 1 : ( 0 | fail);
  derXik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (derXik == NULL) ? 1 : ( 0 | fail);
  tq = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (tq == NULL) ? 1 : ( 0 | fail);
  tx = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (tx == NULL) ? 1 : ( 0 | fail);
  tqp = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (tqp == NULL) ? 1 : ( 0 | fail);
  nQh = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (nQh == NULL) ? 1 : ( 0 | fail);
  dQ = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  fail = (dQ == NULL) ? 1 : ( 0 | fail);


  if (fail)
    return OO_MEMORY;
  /* end - allocate memory */

  /* further initialization of local variables */

  modelica_real diffQ = 0.0, dTnextQ = 0.0, nextQ = 0.0;
  for (i = 0; i < STATES; i++)
  {
    dQ[i] = deltaQFactor * data->modelData->realVarsData[i].attribute.nominal;
    tx[i] = tq[i] = simInfo->startTime;
    qik[i] = state[i];
    xik[i] = state[i];
    derXik[i] = stateDer[i];
    //printf("state[i]: %f\n", state[i]);
    //printf("stateDer[i]: %f\n", stateDer[i]);
    //printf("antonDer[i]: %f\n", antonDer);
    //printf("antonDer[i]: %f\n", antonDer);

    retValue = deltaQ(data, dQ[i], i, &dTnextQ, &nextQ, &diffQ);
    if (OK != retValue)
      return retValue;
    //printf("dTnextQ: %f\n", dTnextQ);
    //printf("nextQ: %f\n", nextQ);
  //  printf("diffQ: %f\n", diffQ);
    tqp[i] = tq[i] + dTnextQ;
    nQh[i] = nextQ;
  }

/* Transform the sparsity pattern into a data structure for an index based access. */
  modelica_integer* der = (modelica_integer*)calloc(ROWS, sizeof(modelica_integer));
  if (NULL==der)
    return OO_MEMORY;
  for (i = 0; i < ROWS; i++)
    der[i] = -1;

  /* how many states are involved in each derivative */
  /* **** This is needed if we have QSS2 or higher **** */
 /* uinteger* numStatesInDer = calloc(ROWS,sizeof(uinteger));
  if (NULL==numStatesInDer) return OO_MEMORY; */
  /*
   *   dx1/dt = x2
   *   dx2/dt = x1 + x2
   *   lead to
   *  StatesInDer[0]-{ 2 }
   *  StatesInDer[1]-{ 1, 2 }
  */
/*  uinteger** StatesInDer = calloc(ROWS,sizeof(uinteger*));
  if (NULL==StatesInDer) return OO_MEMORY; */

  /* count number of states in each derivative*/
/*  for (i = 0; i < pattern->leadindex[ROWS-1]; i++) numStatesInDer[ pattern->index[i] ]++; */
  /* collect memory for all stateindices */
/*  for (i = 0; i < ROWS; i++)
  {
    *(StatesInDer + i) = calloc(numStatesInDer[i], sizeof(uinteger));
    if (NULL==*(StatesInDer + i)) return OO_MEMORY;
  }

  retValue = getStatesInDer(pattern->index, pattern->leadindex, ROWS, STATES, StatesInDer);
  if (OK != retValue) return retValue; */

/*  retValue = getDerWithStateK(pattern->index, pattern->leadindex, der, &numDer, 0);
  if (OK != retValue) return retValue; */
/* End of transformation */

#ifdef D
  FILE* fid=NULL;
  fid = fopen("log_qss.txt","w");
#endif

/* *********************************************************************************** */

/***** Start main simulation loop *****/

  //ek gaan net een of twee iterasies doen om te verstaan wat hier aangaan.
  //int antonInt=0;
  //while(antonInt<4){
  //
  while(solverInfo->currentTime < simInfo->stopTime)
  {
    modelica_boolean syncStep = 0;
    //antonInt++;
    //printf("\tlastdesiredStep: %f\n", solverInfo->lastdesiredStep);
    //printf("\tstateEvents: %d\n", solverInfo->stateEvents);
    //printf("\tdidEventStep: %d\n", solverInfo->didEventStep);
  //  printf("\teventLst.itemSize %d\n", solverInfo->eventLst->itemSize);
  //  printf("\teventLst.length %d\n", solverInfo->eventLst.length);
  //  printf("\tsampleEvents: %d\n", solverInfo->sampleEvents);
    //printf("\tintegratorSteps: %d\n", solverInfo->integratorSteps);
    //antonInt++;
    modelica_integer success = 0;
    threadData->currentErrorStage = ERROR_SIMULATION;
    omc_alloc_interface.collect_a_little();

#if !defined(OMC_EMCC)
    /* try */
    MMC_TRY_INTERNAL(simulationJumpBuffer)
    {
#endif

#ifdef USE_DEBUG_TRACE
    if (useStream[LOG_TRACE])
      printf("TRACE: push loop step=%u, time=%.12g\n", currStepNo, solverInfo->currentTime);
#endif

#ifdef D
    fprintf(fid,"t = %.8f\n",solverInfo->currentTime);
    fprintf(fid,"%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n","tx","x","dx","tq","q","tqp");
    for (i = 0; i < STATES; i++)
    {
      fprintf(fid,"%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\n",tx[i],xik[i],derXik[i],tq[i],qik[i],tqp[i]);
    }
#endif

    currStepNo++;

    ind = minStep(tqp, STATES);

    if (isnan(tqp[ind]))
    {
#ifdef D
      fprintf(fid,"Exit caused by #QNAN!\tind=%d",ind);
#endif
      return ISNAN;
    }
    if (isinf(tqp[ind]))
    {
      /* If all derivatives are zero, the states stay constant and only the
       * time propagates till stop->time.
       */
      warningStreamPrint(LOG_STDOUT, 0, "All derivatives are zero at time %f!.\n", sData->timeValue);
      solverInfo->currentTime = simInfo->stopTime;
      sData->timeValue = solverInfo->currentTime;

      continue;
    }

    qik[ind] = nQh[ind];

    xik[ind] = qik[ind];
    state[ind] = qik[ind];

    tx[ind] = tqp[ind];
    tq[ind] = tqp[ind];

    solverInfo->currentTime = tqp[ind];

#ifdef D
    fprintf(fid,"Index: %d\n\n",ind);
#endif

    /* the state[ind] will change again in dTnextQ*/
    retValue = deltaQ(data, dQ[ind], ind, &dTnextQ, &nextQ, &diffQ);
    if (OK != retValue)
      return retValue;
    tqp[ind] = tq[ind] + dTnextQ;
    nQh[ind] = nextQ;

    data->callback->functionDAE(data, threadData);
    syncStep = simulationUpdate(data, threadData, solverInfo);

    if (0 != strcmp("ia", data->simulationInfo->outputFormat))
    {
      communicateStatus("Running", (solverInfo->currentTime-simInfo->startTime)/(simInfo->stopTime-simInfo->startTime));
    }

    /* get the derivatives depending on state[ind] */
    for (i = 0; i < ROWS; i++)
      der[i] = -1;
    retValue = getDerWithStateK(pattern->index, pattern->leadindex, der, &numDer, ind);


    //data->callback->functionDAE(data, threadData);
    //syncStep = simulationUpdate(data, threadData, solverInfo);


    uinteger k = 0, j = 0;
    for (k = 0; k < numDer; k++)
    {
      j = der[k];
      if (j != ind)
      {
        //printf("xik[j]: %f\n", xik[j]);
        //printf("derXik[j]: %f\n", derXik[j]);
        //printf("solverInfo->currentTime: %f\n", solverInfo->currentTime);
        //printf("tx[j]: %f\n", tx[j]);
        xik[j] = xik[j] + derXik[j] * (solverInfo->currentTime - tx[j]);
        state[j] = xik[j];
        //printf("state[j]: %f\n", state[j]);
        tx[j] = solverInfo->currentTime;
      }
    }
    //syncStep = simulationUpdate(data, threadData, solverInfo);
    //data->callback->functionDAE(data, threadData);
    /*
     * Recalculate all equations which are affected by state[ind].
     * Unfortunately all equations will be calculated up to now. And we need to evaluate
     * the equations as f(t,q) and not f(t,x). So all states were saved onto a local stack
     * and overwritten by q. After evaluating the equations the states are written back.
     */
    for (i = 0; i < STATES; i++)
    {
      xik[i] = state[i];   /* save current state */
      state[i] = qik[i];  /* overwrite current state for dx/dt = f(t,q) */
      //printf("xik[i]: %f\n", xik[i]);
      //printf("state[i]: %f\n", state[i]);
    }
    //data->callback->functionDAE(data, threadData);
    //syncStep = simulationUpdate(data, threadData, solverInfo);

    /* update continous system */
    sData->timeValue = solverInfo->currentTime;
    externalInputUpdate(data);
    data->callback->input_function(data, threadData);
    data->callback->functionODE(data, threadData);
    data->callback->functionAlgebraics(data, threadData);
    data->callback->output_function(data, threadData);
    data->callback->function_storeDelayed(data, threadData);
    data->callback->functionDAE(data, threadData);

    for (i = 0; i < STATES; i++)
    {
      state[i] = xik[i];  /* restore current state */
    }


    /*
     * Get derivatives affected by state[ind] and write back ALL derivatives. After that we have
     * states and derivatives for different times tx.
    */
    //syncStep = simulationUpdate(data, threadData, solverInfo);
    //data->callback->functionDAE(data, threadData);

    for (k = 0; k < numDer; k++)
    {
      j = der[k];
      derXik[j] = stateDer[j];
    }
    derXik[ind] = stateDer[ind];  /* not in every case part of the above derivatives */

    for (i = 0; i < STATES; i++)
    {
      stateDer[i] = derXik[i];  /* write back all derivatives */
    }
    //syncStep = simulationUpdate(data, threadData, solverInfo);
    /* recalculate the time of next change only for the affected states */
    for (k = 0; k < numDer; k++)
    {
       j = der[k];
      retValue = deltaQ(data, dQ[j], j, &dTnextQ, &nextQ, &diffQ);
      if (OK != retValue)
        return retValue;
      tqp[j] = solverInfo->currentTime + dTnextQ;
      nQh[j] = nextQ;
    }
  //syncStep = simulationUpdate(data, threadData, solverInfo);
    /*sData->timeValue = solverInfo->currentTime;*/
    //solverInfo->laststep = solverInfo->currentTime;
    //syncStep = simulationUpdate(data, threadData, solverInfo);
    //printf("syncStep: %d\n", syncStep);
    sim_result.emit(&sim_result, data, threadData);

    //vir events
    //syncStep = simulationUpdate(data, threadData, solverInfo);



    /* check if terminate()=true */
    if (terminationTerminate)
    {
      printInfo(stdout, TermInfo);
      fputc('\n', stdout);
      infoStreamPrint(LOG_STDOUT, 0, "Simulation call terminate() at time %f\nMessage : %s", data->localData[0]->timeValue, TermMsg);
      simInfo->stopTime = solverInfo->currentTime;
    }

    /* terminate for some cases:
     * - integrator fails
     * - non-linear system failed to solve
     * - assert was called
     */
    if (retValIntegrator)
    {
      retValue = -1 + retValIntegrator;
      infoStreamPrint(LOG_STDOUT, 0, "model terminate | Integrator failed. | Simulation terminated at time %g", solverInfo->currentTime);
      break;
    }
    else if (check_nonlinear_solutions(data, 0))
    {
      retValue = -2;
      infoStreamPrint(LOG_STDOUT, 0, "model terminate | non-linear system solver failed. | Simulation terminated at time %g", solverInfo->currentTime);
      break;
    }
    else if (check_linear_solutions(data, 0))
    {
      retValue = -3;
      infoStreamPrint(LOG_STDOUT, 0, "model terminate | linear system solver failed. | Simulation terminated at time %g", solverInfo->currentTime);
      break;
    }
    else if (check_mixed_solutions(data, 0))
    {
      retValue = -4;
      infoStreamPrint(LOG_STDOUT, 0, "model terminate | mixed system solver failed. | Simulation terminated at time %g", solverInfo->currentTime);
      break;
    }
    success = 1;
#if !defined(OMC_EMCC)
    }
    /* catch */
    MMC_CATCH_INTERNAL(simulationJumpBuffer)
#endif
    if (!success)
    {
      retValue =  -1;
      infoStreamPrint(LOG_STDOUT, 0, "model terminate | Simulation terminated by an assert at time: %g", data->localData[0]->timeValue);
      break;
    }

    TRACE_POP /* pop loop */
  }
  /* End of main loop */

#ifdef D
  fprintf(fid,"t = %.8f\n",solverInfo->currentTime);
  fprintf(fid,"%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n","tx","x","dx","tq","q","tqp");
  for (i = 0; i < STATES; i++)
  {
    fprintf(fid,"%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\n",tx[i],xik[i],derXik[i],tq[i],qik[i],tqp[i]);
  }
  fclose(fid);
#endif

  /* free memory*/
   free(der);
 /*  for (i = 0; i < ROWS; i++) free(*(StatesInDer + i));
   free(StatesInDer);
   free(numStatesInDer); */
   free(qik);
   free(xik);
   free(derXik);
   free(tq);
   free(tx);
   free(tqp);
   free(nQh);
   free(dQ);
   /* end - free memory */

  TRACE_POP
  return retValue;
}


/*! static int deltaQ( DATA* data, const modelica_integer index, modelica_real* dTnextQ, modelica_real* nextQ, modelica_real* diffQ)
 *  \brief  Computes the next step in time and quantity for state[index].
 *  \param [ref] [data]  Global data object.
 *  \param [in]  [dQ] Change of quantity for state[index], (nominal value) * 10^-4.
 *  \param [in] [index]  A new step will be computed for state[index].
 *  \param [out] [dTnextQ]  The state will change after dTnextQ second.
 *  \param [out] [nextQ]  Next quantity reached by the state.
 *  \param [out] [diffQ]  Difference between the states current and future value.
 *  \return  [0]  Everything is fine.
 */
static modelica_integer deltaQ( DATA* data, const modelica_real dQ, const modelica_integer index, modelica_real* dTnextQ, modelica_real* nextQ, modelica_real* diffQ)
{

  /* localData[0] because old values in the ringbuffer are not stored in QSS1 and so the ringbuffer will not be rotated. */
  SIMULATION_DATA *sDataOld = (SIMULATION_DATA*)data->localData[0];
  modelica_real* stateDer = sDataOld->realVars + data->modelData->nStates;


  if (stateDer[index] >= 0 )    /* quantity of the state will increase */
  {
    *nextQ = (floor( sDataOld->realVars[index] / dQ ) + 1 ) * dQ;
    if (*nextQ <= (sDataOld->realVars[index] + EPS))
      *nextQ = *nextQ + dQ;
  }
  else
  {
    *nextQ = floor( sDataOld->realVars[index] / dQ ) * dQ;
    if (*nextQ >= (sDataOld->realVars[index] - EPS ))
      *nextQ = *nextQ - dQ;
  }

  *diffQ = fabs(*nextQ - sDataOld->realVars[index]);
  *dTnextQ = fabs(*diffQ / stateDer[index]);

  return OK;
}

/*! static int getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, int* der, unsigned int* numDer, const unsigned int k)
 *  \brief  Returns the indices of all derivatives with state k inside.
 *  \param [ref] [index]
 *  \param [ref] [leadindex]
 *  \param [out] [der]  Derivatives which are influenced by state k.
 *  \param [out] [numDer] Number of influenced derivatives.
 *  \param [in] [k]  State to look for.
 *  \return [0]  Everything is fine.
 */
static modelica_integer getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k)
{
  uinteger start = 0;
  if (0 < k)
    start = leadindex[k - 1];
  uinteger j = 0;
  uinteger i = 0;
  //printf("k: %d\n", k);
  //printf("start: %d\n", start);
  //printf("leadindex[k]: %d\n", leadindex[k]);
  for (i = start; i < leadindex[k]; i++)
  {
    der[j] = index[i];
    //printf("der[j]: %d\n", der[j]);
    j++;
  }
  *numDer = j;
  return OK;
}
/*! static int getStatesInDer(const unsigned int* index, const unsigned int* leadindex, const unsigned int ROWS, const unsigned int STATES, unsigned int** StatesInDer)
 *  \brief  Return the indices of all states in each derivative for an indexed access.
 *  \param [ref] [index]
 *  \param [ref] [leadindex]
 *  \param [in]  [ROWS] number of derivatives
 *  \param [in]  [STATES] number of states
 *  \param [out] [StatesInDer]  index of states in each derivative
 *
 */
static modelica_integer getStatesInDer(const unsigned int* index, const unsigned int* leadindex, const uinteger ROWS, const uinteger STATES, uinteger** StatesInDer)
{
  uinteger i = 0, k = 0; /* loop var */
  uinteger numDer = 0;
  modelica_integer* der = (modelica_integer*)calloc(ROWS, sizeof(modelica_integer));
  uinteger* stackPointer = (uinteger*)calloc(ROWS, sizeof(uinteger));

  if (NULL == der)
    return OO_MEMORY;

  for (i = 0; i < ROWS; i++)
    der[i] = -1;

  for (i = 0; i < ROWS; i++)
    stackPointer[i] = 0;

  /*    Ask for all states in which derivative they occur. */
  for (k = 0; k < STATES; k++)
  {
    getDerWithStateK(index, leadindex, der, &numDer, k);
    for (i = 0; i < ROWS; i++)
    {
      if (der[i] < 0)
        continue;
      /* stackPointer refers to the next free position for der[i] in StatesInDer */
      StatesInDer[ der[i] ][ stackPointer[ der[i] ] ] = k;
      stackPointer[ der[i] ]++;
      der[i] = -1;  // clear all
    }
  }

  free(der);
  free(stackPointer);
  return OK;
}


/*! static unsigned int minStep(const modelica_real* tqp, const unsigned int size )
 *  \brief  Finds the index of the state which will change first.
 *  \param [in] [tqp]  State[i] will change in time tqp[i].
 *  \param [in] [size]  Number of states.
 *  \return  Index of the state which will change first.
 */
static uinteger minStep(const modelica_real* tqp, const uinteger size )
{
  uinteger i = 0;
  uinteger ind = i;
  modelica_real tmin =
#if defined(_MSC_VER)
      NAN;
#else
      1.0/0.0; /* We can have a QNAN at any index and tqp[i] < QNAN will fail in every case. */
#endif

  for (i = 0; i < size; i++)
  {
    if (tqp[i] < tmin && !isnan(tqp[i]))
    {
      ind = i;
      tmin = tqp[i];
    }
  }
  return ind;
}
