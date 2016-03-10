/*
 * By Anton de Villiers
 * HealthQ Technologies
 * Last Updated: 24 February 2016
 * This solver is based on the LIQSS solver proposed by G Migoni & E Kofman in
 * Linearly Implicit Discrete Event Methods for Stiff ODE's
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
enum LIQSS_error_msg
{
	LIQSS_ISNAN = -3L,      /*!< Time of next change is #QNAN. */
	LIQSS_UNKNOWN = -2L,    /*!< Unspecific error. */
	LIQSS_OO_MEMORY = -1L,  /*!< Allocation of memory fails. */
	LIQSS_OK = 0L           /*!< Everything is fine. */
};

const modelica_real LIQSS_EPS = 1e-6;

// Defining boolean in C

#define bool int
#define true 1
#define false 0

// Debugging flag

const bool DEBUG_LIQSS = false;

// Adding an iteration limit

const bool ITERATION_LIMIT = false;
const uinteger ITERATION_LIMIT_VALUE = 150000;
bool LIMIT = false;

// Step-size factor. This is in relation to the nominal value of each state variable.

const modelica_real deltaQFactorLIQSS =0.0001;

static modelica_real calculateQ(const modelica_real der, const modelica_real xik, const modelica_real qLower, const modelica_real qUpper, const modelica_real qChosen);
static modelica_real calculateQLower(const modelica_real qLower, const modelica_real xik, const modelica_real dQ);
static modelica_real nextTime(const modelica_real dQ, const modelica_real der);
static void LIQSS_getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k);
static uinteger minimumStep(const modelica_real* tqp, const uinteger size );
static uinteger calculateState(DATA* data, threadData_t *threadData);

/*
 *  \param [ref] [data]
 *  \param [ref] [solverInfo]
 *
 *  This function performs the simulation controlled by solverInfo.
 */
modelica_integer prefixedName_LIQSSSimulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo){
    TRACE_PUSH

    // This is used for time management
    SIMULATION_INFO *simInfo = data->simulationInfo;
    solverInfo->currentTime = simInfo->startTime;

    uinteger currStepNo = 0, ind = 0;
    modelica_integer retValue = 0;
    modelica_real dTnextQ = 0.0, nextQ = 0.0;

    /*
     * Investigate why this function must be called.
     * Without this, the solver crashes during execution.
     */



    if (data->callback->initialAnalyticJacobianA(data, threadData)){
		infoStreamPrint(LOG_STDOUT, 0, "Jacobian or sparse pattern is not generated or failed to initialize.");
		//return UNKNOWN;
    }

    /* Initialization */
    uinteger i, k, j, eventOccurQ; /* loop var */
    SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
    modelica_real* state = sData->realVars;
    modelica_real* stateDer = sData->realVars + data->modelData->nStates;
    const SPARSE_PATTERN* pattern = &(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern);
    const uinteger ROWS = data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sizeRows;
    const uinteger STATES = data->modelData->nStates;
    uinteger numDer = 0;  /* number of derivatives influenced by state k */
	modelica_boolean fail = 0;  /* When the solver fails */
	modelica_real* qChosen = NULL;  /* Approximation of states */
	modelica_real* qLower = NULL;  /* q__ j (t) */
	modelica_real* qUpper = NULL;  /* q^-- j (t) */
	modelica_real* xik = NULL;  /* state variable values */
	modelica_real* time = NULL;    /* Time of approximation, because not all approximations are calculated at a specific time, each approx. has its own timestamp */
	modelica_real* timeOld = NULL;  /* The last time the corresponding state variable was updated. */
	modelica_real* dQ = NULL;    /* change in quantity of every state, default = deltaQFactor */
	modelica_integer* der = NULL;

    /* allocate memory*/
    qChosen = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (qChosen == NULL) ? 1 : ( 0 | fail);
  	xik = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (xik == NULL) ? 1 : ( 0 | fail);
  	qLower = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (qLower == NULL) ? 1 : ( 0 | fail);
  	qUpper = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (qUpper == NULL) ? 1 : ( 0 | fail);
  	time = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (time == NULL) ? 1 : ( 0 | fail);
  	timeOld = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (timeOld == NULL) ? 1 : ( 0 | fail);
  	dQ = (modelica_real*)malloc(STATES * sizeof(modelica_real));
  	fail = (dQ == NULL) ? 1 : ( 0 | fail);
    der = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
    fail = (der == NULL) ? 1 : ( 0 | fail);

    if (fail)
        return OO_MEMORY;
    /* end - allocate memory */

    for (i = 0; i < STATES; i++){
		dQ[i] = deltaQFactorLIQSS * data->modelData->realVarsData[i].attribute.nominal;
		//dQ[i] = deltaQFactorLIQSS;
    	time[i] = timeOld[i] = simInfo->startTime;
		xik[i] = state[i];
		qLower[i] = state[i] - dQ[i];
		qUpper[i] = state[i] + dQ[i];

		if(stateDer[i] > 0)
			qChosen[i] = qUpper[i];
		else if(stateDer[i] < 0)
			qChosen[i] = qLower[i];
		else
			qChosen[i] =state[i];
		dTnextQ = nextTime(dQ[i], stateDer[i]);
		time[i] = time[i] + dTnextQ;
		if(DEBUG_LIQSS){
			printf("%f\t%d\tdQ[i]: %.12f\n", solverInfo->currentTime, i, dQ[i]);
			printf("%f\t%d\txik[i]: %.10f\n", solverInfo->currentTime, i, xik[i]);
			printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
			printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
			printf("%f\t%d\tdTnextQ: %.12f \n", solverInfo->currentTime, i, time[i]);
			printf("%f\t%d\tqChosen[i]: %.12f \n", solverInfo->currentTime, i, qChosen[i]);
		}
    }
    currStepNo++;


	/* get the derivatives depending on state[ind] */
	for (i = 0; i < ROWS; i++)
		der[i] = -1;
	for (k = 0; k < STATES; k++){
		ind =k;
		LIQSS_getDerWithStateK(pattern->index, pattern->leadindex, der, &numDer, ind);
		printf("numDer %d \n", numDer);
		printf("ind %d \n", ind);
		for(i =0;i<ROWS;i++){
			printf("der %d: %d\n",i,der[i]);
		}
	}

	printf("What is difference between STATES %d and ROWS %d\n\n", STATES, ROWS);

	i=0;
	int start =0;
	for(k =0;k<ROWS;k++){
		printf("ROW: %d\n",k);
		if(k!=0){
			start=pattern->leadindex[k-1];
			printf("** %d\n", pattern->leadindex[k-1]);
			printf("** %d\n", pattern->leadindex[k]);
		}
		for (i=start;i<pattern->leadindex[k];i++){
			printf("-- %d\n",pattern->index[i]);
		}
	}


//	static void LIQSS_getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k){
//		uinteger start = 0, j = 0;
//		if (0 < k)
//			start = leadindex[k-1];
//		for (uinteger i = start; i < leadindex[k]; i++){
//			der[j] = index[i];
//			j++;
//		}
//		*numDer = j;
//	}

//	for(k=0;k<(sizeof(pattern->index) / sizeof(int));k++){
//		printf("gfgfg\n");
//	}



    /* Find the next time step before going into the loop */
    ind = minimumStep(time, STATES);


    /***** Start main simulation loop *****/
    while(solverInfo->currentTime < simInfo->stopTime && !LIMIT){
    	if(ITERATION_LIMIT)
    		if(currStepNo == ITERATION_LIMIT_VALUE - 1)
    			LIMIT = true;

		if(DEBUG_LIQSS){
			printf("\\***************************************************\\\n");
		}
		currStepNo++;

		//Capping the time increase to dQ[ind]
		if(time[ind]-solverInfo->currentTime>dQ[ind])
			time[ind]=dQ[ind]+ solverInfo->currentTime;

		if(time[ind] <= simInfo->stopTime)
			solverInfo->currentTime = time[ind];
		else
			solverInfo->currentTime = simInfo->stopTime;

		for (k = 0; k < STATES; k++){
		    //j = der[k];
			j=k;
			xik[j] = xik[j] + stateDer[j] * (solverInfo->currentTime - timeOld[j]);
			timeOld[j] = solverInfo->currentTime;
			if(DEBUG_LIQSS){
				printf("%f %d\txik[j]: %.12f \n", solverInfo->currentTime, j, xik[j]);
				printf("%f %d\tstateDer[j]: %.12f \n", solverInfo->currentTime, j, stateDer[j]);
				printf("%f %d\ttime[j]: %.12f \n", solverInfo->currentTime, j, time[j]);
				printf("%f %d\tqChosen[j]: %.12f \n", solverInfo->currentTime, j, qChosen[j]);
					if(i==STATES-1)
						printf("\n");
			}
		}


//		/*
//		 * We have to update here since we will use the xik values,
//		 * which are set as the state variables when updating qChosen
//		 * for all the required state variables.
//		 */
//
		for(i =0;i<STATES;i++){
			state[i] =qChosen[i];
		}

		eventOccurQ=calculateState(data, threadData);
        simulationUpdate(data, threadData, solverInfo);

		//numDer
		for (j = 0; j < STATES; j++){
    	    //k = der[j];
			k=j;
			/*
			* Only alter the q values for the state variable whose q values has not yet reached the x values
			*/
			if(1){
				qLower[k] = calculateQLower(qLower[k], xik[k], dQ[k]);
				qUpper[k] = qLower[k] + 2*dQ[k];

				/*
				* To perform calculateQ you may have to calculate qTilde.
				* In this event, we will have to have to know fqUpperJ and fqLowerJ
				* We won't necessarily require these values, but I am calculating them either way.
				* Can tidy this up perhaps -- let's do it now...
				*/

				state[k] = qUpper[k];
				calculateState(data, threadData);
				modelica_real fqUpperJ = stateDer[k];

				state[k] = qLower[k];
				calculateState(data, threadData);
				modelica_real fqLowerJ = stateDer[k];

				state[k] = qChosen[k];
				calculateState(data, threadData);

				if ((fqLowerJ*(qLower[k] - xik[k])) >= 0)
					qChosen[k] = qLower[k];
				else if ((fqUpperJ*(qUpper[k] - xik[k])) >= 0 && (fqLowerJ*(qLower[k] - xik[k])) < 0)
					qChosen[k] = qUpper[k];
				else{
					modelica_real Ajj = (fqUpperJ - fqLowerJ)/(qUpper[k] - qLower[k]);
					if(Ajj == 0)
						qChosen[k] = qChosen[k];
					else
					    qChosen[k] = qUpper[k] - (fqUpperJ/Ajj);
				}


				//if(DEBUG_LIQSS){
				if(DEBUG_LIQSS || qChosen[k] < qLower[k] || qChosen[k] > qUpper[k]){

					printf("%f\t%d\tfqLowerJ: %.12f \n", solverInfo->currentTime, k, fqLowerJ);
					printf("%f\t%d\tfqUpperJ: %.12f \n", solverInfo->currentTime, k, fqUpperJ);
					printf("%f\t%d\tqLower[i]: %.12f \n", solverInfo->currentTime, k, qLower[k]);
					printf("%f\t%d\tqUpper[i]: %.12f \n", solverInfo->currentTime, k, qUpper[k]);
					printf("%f\t%d\tqChosen[i]: %.12f \n", solverInfo->currentTime, k, qChosen[k]);
				}
			}
			else{
				if(DEBUG_LIQSS)
					printf("Indeks waarvoor qChosen[k] nie herbereken word nie is: %d\n",k);
			}
		}

        /* No we must calculate the derivatives but with the qChosen values.
         * Another fat cheat is used here. We update the global system with
         * the qChosen values. We will rewrite the actual global variables later.
         * */

//		for (i = 0; i < STATES; i++){
//			//k = der[i];
//			k=i;
//			state[k] = qChosen[k];
//			if(DEBUG_LIQSS){
//				if(i==STATES-1)
//					printf("\n");
//				printf("%f\t%d\tqChosen[i]: %.12f \n", solverInfo->currentTime, k, qChosen[k]);
//				if(i==STATES-1)
//					printf("\n");
//			}
//		}

        /* update continuous system */
//        sData->timeValue = solverInfo->currentTime;
//        calculateState(data, threadData);
//        simulationUpdate(data, threadData, solverInfo);


         /* Update the global variables for this round
          * Remember we have calculated qChosen earlier
          * These values has to be captured for the current time step
          * */

         for (i = 0; i < STATES; i++){
             state[i] = qChosen[i];
         }





         /* update continuous system */
         sData->timeValue = solverInfo->currentTime;
         eventOccurQ=calculateState(data, threadData);
         simulationUpdate(data, threadData, solverInfo);
         sim_result.emit(&sim_result, data, threadData);






         /* Now we calculate the next time */
		for (i = 0; i < STATES; i++){


			dTnextQ = nextTime(dQ[i], stateDer[i]);
			time[i] = sData->timeValue + dTnextQ;
			if(DEBUG_LIQSS){
				if(i==0)
					printf("We have now assigned qChosen to state. We did update. Now calculating nextTime.\n");
				printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
				printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
				printf("%f\t%d\txik[i]: %.10f\n", solverInfo->currentTime, i, xik[i]);
				printf("%f\t%d\tqChosen[i]: %.10f\n", solverInfo->currentTime, i, qChosen[i]);
				printf("%f\t%d\tdTnextQ[i]: %.10f\n", solverInfo->currentTime, i, dTnextQ);
				printf("%f\t%d\ttime[i]: %.10f\n", solverInfo->currentTime, i, time[i]);
				if(i==STATES-1)
					printf("\n");
			}
		}

		for (i = 0; i < STATES; i++){
			state[i] = qChosen[i];
		}



		sData->timeValue = solverInfo->currentTime;
		eventOccurQ=calculateState(data, threadData);
		simulationUpdate(data, threadData, solverInfo);



		/*
		* In the case where a state event is triggered, the state variables are resetted.
		* We also have to reset the q values for the reinitialisation process.
		*/

		if(eventOccurQ){
			for (i = 0; i < STATES; i++){
				timeOld[i] = solverInfo->currentTime;
				xik[i] = state[i];
				qLower[i] = state[i] - dQ[i];
				qUpper[i] = state[i] + dQ[i];
				if(stateDer[i] > 0)
					qChosen[i] = qUpper[i];
				else if(stateDer[i] < 0)
					qChosen[i] = qLower[i];
				else
					qChosen[i] =state[i];
				if(DEBUG_LIQSS){
					printf("When an event occurs\n\n");
					printf("%f\t%d\tstateDer[i]: %.10f\n", solverInfo->currentTime, i, stateDer[i]);
					printf("%f\t%d\tstate[i]: %.10f\n", solverInfo->currentTime, i, state[i]);
					printf("%f\t%d\txik[i]: %.10f\n", solverInfo->currentTime, i, xik[i]);
					printf("%f\t%d\tqChosen[i]: %.10f\n", solverInfo->currentTime, i, qChosen[i]);
					printf("%f\t%d\tdTnextQ[i]: %.10f\n", solverInfo->currentTime, i, dTnextQ);
					printf("%f\t%d\ttime[i]: %.10f\n", solverInfo->currentTime, i, time[i]);
				}
			}
		}

		ind = minimumStep(time, STATES);



		sim_result.emit(&sim_result, data, threadData);

		TRACE_POP /* pop loop */
    }
    /* End of main loop */

    /* free memory*/
	free(qChosen);
	free(qLower);
	free(xik);
	free(qUpper);
	free(time);
	free(timeOld);
	free(dQ);
	free(der);
	/* end - free memory */

	printf("The number of iterations were: %d\n", currStepNo);

	TRACE_POP
	return retValue;
}


/*
 *  \brief  Returns the indices of all derivatives with state k inside.
 *  \param [ref] [index]
 *  \param [ref] [leadindex]
 *  \param [out] [der]  Derivatives which are influenced by state k.
 *  \param [out] [numDer] Number of influenced derivatives.
 *  \param [in] [k]  State to look for.
 */

static void LIQSS_getDerWithStateK(const unsigned int *index, const unsigned int* leadindex, modelica_integer* der, uinteger* numDer, const uinteger k){
	uinteger start = 0, j = 0;
	if (0 < k)
		start = leadindex[k-1];
	for (uinteger i = start; i < leadindex[k]; i++){
		der[j] = index[i];
		j++;
	}
	*numDer = j;
}

/*
 *  \brief  Finds the index of the state which will change first.
 *  \param [in] [tqp]  State[i] will change in time tqp[i].
 *  \param [in] [size]  Number of states.
 *  \return  Index of the state which will change first.
 */

static uinteger minimumStep(const modelica_real* tqp, const uinteger size){
	uinteger ind = 0;
	modelica_real tmin =
	#if defined(_MSC_VER)
		NAN;
	#else
		1.0/0.0;
	#endif
	for (uinteger i = 0; i < size; i++){
		if (tqp[i] < tmin && !isnan(tqp[i])) {
			ind = i;
		tmin = tqp[i];
		}
	}
	return ind;
}


/*
 *  \brief  Calculates the next time period for a specific state variable
 *  \param [in] [dQ]  deltaQ for state variable.
 *  \param [in] [der]  derivative for state variable.
 *  \return  The time when x will reach the chosen q value.
 */

static modelica_real nextTime(const modelica_real dQ, const modelica_real der){
	if (der == 0)
		return dQ/1e-2;
	else
		return dQ/fabs(der);
}

// Calculates the qLower value for a specific state variable

/*
 *  \brief  Calculates the qLower value of a state variable
 *  \param [in] [qLower]  the previous qLower value.
 *  \param [in] [xik]  value of the state variable.
 *  \param [in] [dQ]  deltaQ for state variable.
 *  \return  The qLower value based on the current time.
 */

static modelica_real calculateQLower(const modelica_real qLower, const modelica_real xik, const modelica_real dQ){
	if(DEBUG_LIQSS){
		printf("Inside calculateQLower\nqLower %f\txik %f\tdq %f\n",qLower,xik,dQ);
		printf("Inside calculateQLower\nxik - qLower %f\t2dQ %f\n",xik - qLower,2*dQ);
	}
	if(xik - qLower <= LIQSS_EPS){
		//printf("This one 1\n");
		return qLower - dQ;

	}
	else if(xik - qLower >= 2*dQ-LIQSS_EPS){
		//printf("This one 2\n");
		return qLower + dQ;

	}
	else{
		//printf("This one 3\n");
		return qLower;

	}
}

/*
 *  \brief  Computes the ODE's and updates the results
 *  \param [in] [data]
 *  \param [in] [threadData]
 */

static uinteger calculateState(DATA* data, threadData_t *threadData){
	externalInputUpdate(data);
	data->callback->input_function(data, threadData);
	data->callback->functionODE(data, threadData);
	data->callback->functionAlgebraics(data, threadData);
	data->callback->output_function(data, threadData);
	data->callback->function_storeDelayed(data, threadData);;
	data->callback->functionDAE(data, threadData);
	if(data->simulationInfo->needToIterate){
		if(DEBUG_LIQSS){
			printf("needToIterate is een\n");
			printf("nTI:  %d\n", data->simulationInfo->needToIterate);
		}
		return 1;
	}
	else
		return 0;
}

