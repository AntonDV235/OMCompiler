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
const modelica_real deltaQFactor =1;

#define bool int
#define true 1
#define false 0

const bool DEBUG_LIQSS = true;


/* Needed if we want to write all the variables into a file*/
/* #define D */

static modelica_real calculateQ(const modelica_real derXik, const modelica_real xik, const modelica_real qLower, const modelica_real qUpper, const modelica_real qChosen);
static modelica_real calculateQLower(const modelica_real qLower, const modelica_real xik, const modelica_real dQ);
static modelica_real nextTime(const modelica_real dQ, const modelica_real derXik);
static modelica_integer deltaQ( DATA* data,const modelica_real dQ, const modelica_integer index, modelica_real* dTnextQ, modelica_real* nextQ, modelica_real* diffQ);
static modelica_integer getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k);
static modelica_integer getStatesInDer(const unsigned int* index, const unsigned int* leadindex, const uinteger ROWS, const uinteger STATES, uinteger** StatesInDer);
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

  uinteger currStepNo = 0;
  modelica_integer retValue = 0;
  uinteger ind = 0;

  solverInfo->currentTime = simInfo->startTime;

  // Hierdie shit is nodig -- crash die solver as dit nie daar is nie.
  // Fok weet hoekom, though. Dit print tog net logs, sekerlik.
  // Miskien is dit dat infoStreamPrint "moet" gebeur in leeftyd van solver.
  // Sal moet investigate.
  if (data->callback->initialAnalyticJacobianA(data, threadData))
  {
    infoStreamPrint(LOG_STDOUT, 0, "Jacobian or sparse pattern is not generated or failed to initialize.");
    //return UNKNOWN;
  }
  //printSparseStructure(data, LOG_SOLVER);

/* *********************************************************************************** */
  /* Initialization */
  uinteger i, k, j = 0; /* loop var */
  SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
  modelica_real* state = sData->realVars;
  modelica_real* stateDer = sData->realVars + data->modelData->nStates;




  const SPARSE_PATTERN* pattern = &(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern);
  const uinteger ROWS = data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sizeRows;
  const uinteger STATES = data->modelData->nStates;
  uinteger numDer = 0;  /* number of derivatives influenced by state k */



  // must delete all of this
    modelica_boolean fail = 0;
    modelica_real* qik = NULL;  /* Approximation of states */
    modelica_real* tq = NULL;    /* Time of approximation, because not all approximations are calculated at a specific time, each approx. has its own timestamp */
    modelica_real* tx = NULL;    /* Time of the states, because not all states are calculated at a specific time, each state has its own timestamp */
    modelica_real* tqp = NULL;    /* Time of the next change in state */
    modelica_real* nQh = NULL;    /* next value of the state */

    modelica_real* qChosen = NULL;  /* Approximation of states */
    modelica_real* qLower = NULL;  /* q__ j (t) */
    modelica_real* qUpper = NULL;  /* q^-- j (t) */
    modelica_real* xik = NULL;  /* states */
    modelica_real* derXik = NULL;  /* Derivative of states */
    modelica_real* time = NULL;    /* Time of approximation, because not all approximations are calculated at a specific time, each approx. has its own timestamp */
    modelica_real* timeOld = NULL;
    modelica_real* dQ = NULL;    /* change in quantity of every state, default = deltaQFactor */
    modelica_integer* der = NULL;


  /* allocate memory*/
    qik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
    fail = (qik == NULL) ? 1 : ( 0 | fail);
    tq = (modelica_real*)calloc(STATES, sizeof(modelica_real));
    fail = (tq == NULL) ? 1 : ( 0 | fail);
    tx = (modelica_real*)calloc(STATES, sizeof(modelica_real));
    fail = (tx == NULL) ? 1 : ( 0 | fail);
    tqp = (modelica_real*)calloc(STATES, sizeof(modelica_real));
    fail = (tqp == NULL) ? 1 : ( 0 | fail);
    nQh = (modelica_real*)calloc(STATES, sizeof(modelica_real));
    fail = (nQh == NULL) ? 1 : ( 0 | fail);

    qChosen = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (qChosen == NULL) ? 1 : ( 0 | fail);
  	xik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (xik == NULL) ? 1 : ( 0 | fail);
  	qLower = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (qLower == NULL) ? 1 : ( 0 | fail);
  	qUpper = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (qUpper == NULL) ? 1 : ( 0 | fail);
  	derXik = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (derXik == NULL) ? 1 : ( 0 | fail);
  	time = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (time == NULL) ? 1 : ( 0 | fail);
  	timeOld = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (timeOld == NULL) ? 1 : ( 0 | fail);
  	dQ = (modelica_real*)calloc(STATES, sizeof(modelica_real));
  	fail = (dQ == NULL) ? 1 : ( 0 | fail);
  der = (modelica_integer*)calloc(ROWS, sizeof(modelica_integer));
  fail = (der == NULL) ? 1 : ( 0 | fail);

  if (fail)
    return OO_MEMORY;
  /* end - allocate memory */

  /* further initialization of local variables */

  modelica_real diffQ = 0.0, dTnextQ = 0.0, nextQ = 0.0;
  for (i = 0; i < STATES; i++){
    dQ[i] = deltaQFactor * data->modelData->realVarsData[i].attribute.nominal;
    time[i] = timeOld[i] = simInfo->startTime;
    qChosen[i] = state[i];
    xik[i] = state[i];
    derXik[i] = stateDer[i];
    qLower[i] = qChosen[i] - dQ[i];
    qUpper[i] = qChosen[i] + dQ[i];
    dTnextQ = nextTime(dQ[i], stateDer[i]);
    time[i] = time[i] + dTnextQ;
    // For dubugging purposes
    if(DEBUG_LIQSS){
		printf("%f\t%d\tdQ[i]: %.12f\n", solverInfo->currentTime, i, dQ[i]);
		printf("%f\t%d\txik[i]: %.10f\n", solverInfo->currentTime, i, xik[i]);
		printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
		printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
    	printf("%f\t%d\tdTnextQ: %.12f \n", solverInfo->currentTime, i, time[i]);
    	printf("%f\t%d\tqChosen: %.12f \n\n", solverInfo->currentTime, i, qChosen[i]);
    }
  }

  currStepNo++;

  /* Find the next time step before going into the loop */
  ind = minStep(time, STATES);


/* *********************************************************************************** */
/***** Start main simulation loop *****/


  while(solverInfo->currentTime < simInfo->stopTime && currStepNo <5)
  {

	if(DEBUG_LIQSS){
		printf("\\*****************************************\\\n");
	}
    modelica_boolean syncStep = 0;

    currStepNo++;

    solverInfo->currentTime = time[ind];



    /* get the derivatives depending on state[ind] */
    for (i = 0; i < ROWS; i++)
      der[i] = -1;
    getDerWithStateK(pattern->index, pattern->leadindex, der, &numDer, ind);




    for (k = 0; k < numDer; k++){
          j = der[k];
    	  xik[j] = xik[j] + stateDer[j] * (solverInfo->currentTime - timeOld[j]);
    	  timeOld[j] = time[j];
    	  time[j] = solverInfo->currentTime;
        }

    for (k = 0; k < STATES; k++){
          qLower[k] = calculateQLower(qLower[k], xik[k], dQ[k]);
          qUpper[k] = qLower[k] + 2*dQ[k];

          /*
           * To perform calculateQ you may have to calculate qTilde.
           * In this event, we will have to have to know fqUpperJ and fqLowerJ
           * We won't necessarily require these values, but I am calculating them either way.
           * Can tidy this up perhaps -- let's do it now...
           *
           * Remove or alter calculateQ later
           */



		  	if ((stateDer[k]*(qLower[k] - xik[k])) >= 0)
		  		qChosen[k] = qLower[k];
		  	else if ((stateDer[k]*(qUpper[k] - xik[k])) >= 0 && (stateDer[k]*(qLower[k] - xik[k])) < 0)
		  		qChosen[k] = qUpper[k];
		  	else{
		  		printf ("Moet nou qTilde bereken!\n");

		  		 state[k] = qUpper[k];
				  externalInputUpdate(data);
				  data->callback->input_function(data, threadData);
				  data->callback->functionODE(data, threadData);
				  data->callback->functionAlgebraics(data, threadData);
				  data->callback->output_function(data, threadData);
				  data->callback->function_storeDelayed(data, threadData);
				  data->callback->functionDAE(data, threadData);
				  modelica_real fqUpperJ = stateDer[k];

				  state[k] = qLower[k];
				  externalInputUpdate(data);
				  data->callback->input_function(data, threadData);
				  data->callback->functionODE(data, threadData);
				  data->callback->functionAlgebraics(data, threadData);
				  data->callback->output_function(data, threadData);
				  data->callback->function_storeDelayed(data, threadData);
				  data->callback->functionDAE(data, threadData);
				  modelica_real fqLowerJ = stateDer[k];

				  modelica_real Ajj = (fqUpperJ - fqLowerJ)/(qUpper[k] - qLower[k]);
				  if(Ajj == 0)
					  qChosen[k] = qChosen[k]; //redundant, maar ja, net om te bevestig vir volledigheidshalwe -- doen nie skade nie
				  else
					  qChosen[k] = qUpper[k] - (fqUpperJ/Ajj);


		  	}



          //qChosen[k] = calculateQ(derXik[k], xik[k], qLower[k], qUpper[k], qChosen[k]);
        }

        /* No we must calculate the derivatives but with the qChosen values.
         * Another fat cheat is used here. We update the global system with
         * the qChosen values. We will rewrite the actual global variables later.
         * */

        for (i = 0; i < STATES; i++){
          state[i] = qChosen[i];
        }

        /* update continuous system */
        sData->timeValue = solverInfo->currentTime;
        externalInputUpdate(data);
        data->callback->input_function(data, threadData);
        data->callback->functionODE(data, threadData);
        data->callback->functionAlgebraics(data, threadData);
        data->callback->output_function(data, threadData);
        data->callback->function_storeDelayed(data, threadData);
        data->callback->functionDAE(data, threadData);
        syncStep = simulationUpdate(data, threadData, solverInfo);
        if(DEBUG_LIQSS){
        	printf("We have now assigned qChosen to state. We did update. Now calculating nextTime.\n");
        	for (i = 0; i < STATES; i++){
				printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
				printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
			}
        	printf("\n");
        }

        /* Now we calculate the next time */

            for (i = 0; i < STATES; i++){
              dTnextQ = nextTime(dQ[i], stateDer[i]);
              //printf("dTnextQ: %f %d\n", dTnextQ, i);
              time[i] = time[i] + dTnextQ;
            }

            ind = minStep(time, STATES);

            /* Update the global variables for this round
             * Remember we have calculated xik earlier
             * These values has to be captured for the current time step
             * */

            // Update state[i] and update the global results

            for (i = 0; i < STATES; i++){
              state[i] = xik[i];
            }

            /* update continuous system */
            sData->timeValue = solverInfo->currentTime;
            externalInputUpdate(data);
            data->callback->input_function(data, threadData);
            data->callback->functionODE(data, threadData);
            data->callback->functionAlgebraics(data, threadData);
            data->callback->output_function(data, threadData);
            data->callback->function_storeDelayed(data, threadData);
            data->callback->functionDAE(data, threadData);
            syncStep = simulationUpdate(data, threadData, solverInfo);


            if(DEBUG_LIQSS){
            	for (i = 0; i < STATES; i++){
					printf("%f\t%d\tdQ[i]: %.12f\n", solverInfo->currentTime, i, dQ[i]);
					printf("%f\t%d\txik[i]: %.10f\n", solverInfo->currentTime, i, xik[i]);
					printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
					printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
					printf("%f\t%d\ttime[i]: %.12f \n", solverInfo->currentTime, i, time[i]);
					printf("%f\t%d\tqChosen: %.12f \n\n", solverInfo->currentTime, i, qChosen[i]);
            	}
            }

    sim_result.emit(&sim_result, data, threadData);


    TRACE_POP /* pop loop */
  }
  /* End of main loop */

  /* free memory*/
   free(der);
   free(qik);
   free(xik);
   free(tq);
   free(tx);
   free(tqp);
   free(nQh);
   free(dQ);
   /* end - free memory */

  TRACE_POP
  return retValue;
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

// Calculates the next time period for a specific state variable

static modelica_real nextTime(const modelica_real dQ, const modelica_real derXik){
	modelica_real der = derXik;
    if (der == 0)
	  der = 1e-5;
	return dQ/fabs(der);
}

// Calculates the qLower value for a specific state variable

static modelica_real calculateQLower(const modelica_real qLower, const modelica_real xik, const modelica_real dQ){
	if(xik - qLower <= 0)
		return qLower - dQ;
	else if(xik - qLower >= 2*dQ)
		return qLower + dQ;
	else
		return qLower;
}

// Calculates the actual q variable for a state variable --  this variable (the selected q value) is named qChosen

static modelica_real calculateQ(const modelica_real derXik, const modelica_real xik,
		const modelica_real qLower, const modelica_real qUpper, const modelica_real qChosen){
	/* Gaan nie qTilde bereken nie -- sal dit later bysit indien nodig  */
	if ((derXik*(qLower - xik)) >= 0)
		return qLower;
	else if ((derXik*(qUpper - xik)) >= 0 && (derXik*(qLower - xik)) < 0)
		return qUpper;
	else{
		printf ("Moet nou qTilde bereken");
		return qChosen;
	}
}
