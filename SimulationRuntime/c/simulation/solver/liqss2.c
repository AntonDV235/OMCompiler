/*
 * By Anton de Villiers
 * HealthQ Technologies
 * Last Updated: 24 February 2016
 * This solver is based on the LIQSS2 solver
 * incorporated from QSS Solver (applicaiton)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solver_main.h"
#include "simulation/simulation_runtime.h"
#include "simulation/results/simulation_result.h"
#include "openmodelica_func.h"
#include "util/omc_error.h"
#include "simulation/options.h"
#include "simulation/solver/liqss2Operations.h"
#include <unistd.h>


int prefixedName_LIQSS2Simulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo){
	TRACE_PUSH
	data->callback->initialAnalyticJacobianA(data, threadData);
    const uinteger STATES = data->modelData->nStates; //Number of state variables
	modelica_real ft = data->simulationInfo->stopTime;
	int iterations = 0;
	modelica_real t = 0, dt =0, elapsed, mpr1 = 0, mpr2 = 0;
	const SPARSE_PATTERN* pattern = &(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern);
	int dt_index = 0;
	uinteger start =0,jIndex = 0, i = 0, j = 0;

	const modelica_real tolerance = 0.000001;
	const modelica_real absTolerance = 0.000000001;

	LIQSS2_quantizerState quantizer;
	quantizer = LIQSS2_QuantizerState ();
	LIQSS2_init(quantizer, STATES);

	for(i=0; i<STATES; i++){
		if(i>0)
			start = pattern->leadindex[i-1];
		quantizer->nSD[i] = pattern->leadindex[i] - start;
		jIndex = 0;
		for(uinteger j = start; j < pattern->leadindex[i]; j++){
			quantizer->SD[i][jIndex] = pattern->index[j];
			jIndex++;
		}
	}

	// Initial values and the initial derivative values
	SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
	modelica_real* state = sData->realVars;
    modelica_real* stateDer = sData->realVars + data->modelData->nStates;
	modelica_real QEvent = false;

	if(LIQSS_write(quantizer, data, 1, t, STATES, dt_index) == 1){
		return 1;
	}

	//LIQSS2_printModelData(data);

	LIQSS_initialiation(quantizer, state, stateDer, STATES, ft, tolerance, absTolerance);

	if(LIQSS_write(quantizer, data, 2, t, STATES, dt_index) == 1){
		return 1;
	}


	printf("\nThe Initialization stage\n");
	dt = quantizer->mpr[0];
	dt_index = 0;
	for(i=0;i<STATES;i++){
		printf("nTime[i]:  %d  %.16f\n", i, quantizer->mpr[i]);
		if(dt > quantizer->mpr[i]){
			dt = quantizer->mpr[i];
			dt_index = i;
		}
	}
	t=dt;
	solverInfo->currentTime = t;


	// Now we do the integration iteratively
	// This is based on the QSS_SEQH_integrate function in qss_seqh_integrator.c
	while (t < ft){
		iterations += 1;
		printf("\nBeginning of iteration\n");
		printf("t: %.16f\n", t);
		printf("dt_index: %d\n", dt_index);
		// printf("t - tx[dt_index]: %f\n", t - tx[dt_index]);
		//printf("x[dt_index][0]: %f\n", x[dt_index][0]);
		//printf("x[dt_index][1]: %f\n", x[dt_index][1]);
		//printf("x[dt_index][2]: %f\n", x[dt_index][2]);
		elapsed = t - quantizer->tx[dt_index];
		printf("Elapsed: %.16f\n", elapsed);
		// advanceTime in src/common/utils.c
		quantizer->x[dt_index][0] = calcState(quantizer->x, dt_index, elapsed);
		quantizer->x[dt_index][1] = calcDerivative(quantizer->x, dt_index, elapsed);
		quantizer->tx[dt_index] = t;
		quantizer->lqu[dt_index] = fabs(quantizer->x[dt_index][0]) * tolerance;
		if(quantizer->lqu[dt_index] < absTolerance){
			quantizer->lqu[dt_index] = absTolerance;
		}

		// QA_updateQuantizedState -- liqss2.c
		LIQSS2_updateQuantizedState(quantizer, t, dt_index);

		// printf("Na QA_updateQuantizedState\n");
		// for(i = 0;i<STATES;i++){
		// 	for(j =0;j<order+1;j++){
		// 		printf("q[i][j]: %d %d   %.16f\n",i,j,q[i][j]);
		// 	}
		// }
		// printf("\n");
		// for(i = 0;i<STATES;i++){
		// 	for(j =0;j<order+1;j++){
		// 		printf("x[i][j]: %d %d   %.16f\n",i,j,x[i][j]);
		// 	}
		// }
		// printf("\n");

		quantizer->tq[dt_index] = t;

		// QA_nextTime -- liqss2.c
		if(quantizer->x[dt_index][2] == 0)
			quantizer->nTime[dt_index] = 1.0 / 0.0;
		else
			quantizer->nTime[dt_index]  = t + sqrt(fabs(quantizer->lqu[dt_index] / quantizer->x[dt_index][2]));

		for(i=0;i<quantizer->nSD[dt_index];i++){

			j= quantizer->SD[dt_index][i];
			elapsed = t - quantizer->tx[j];
			if(elapsed > 0){
				quantizer->x[j][0] = calcState(quantizer->x,j,elapsed);
				quantizer->tx[j] = t;
			}

			// FRW_recomputeDerivatives -- SYM_recomputeDerivatives qss_frw_imp.c
			elapsed = t - quantizer->tq[j];
			if(elapsed > 0){
				quantizer->q[j][0] = calcQState(quantizer->q,j,elapsed);
				quantizer->tq[j] = t;
			}
		}

		// FRW_recomputeDerivatives -- SYM_recomputeDerivatives qss_frw_imp.c

		SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
	    modelica_real* state = sData->realVars;
	    modelica_real* stateDer = sData->realVars + data->modelData->nStates;
		uinteger i = 0;
		for(i = 0; i< data->modelData->nStates; i++){
			state[i] = quantizer->q[i][0];
		}
		sData->timeValue = solverInfo->currentTime;
		QEvent = LIQSS2_calculateState(data, threadData);
		simulationUpdate(data, threadData, solverInfo);
		sim_result.emit(&sim_result, data, threadData);
		for(i = 0; i< data->modelData->nStates; i++){
			quantizer->x[i][1] = stateDer[i];
		}
		for(i = 0; i< STATES; i++){
			quantizer->x[i][2] = ddx(quantizer->q,i,state);
		}

		// FRW_nextEventTime -- still needs to investigate

		if(QEvent){
			for(i=0;i<STATES;i++){
				quantizer->tx[i] = t;
				quantizer->tq[i] = t;
				quantizer->ltq[i] = t;
				quantizer->lt[i] = t;
				quantizer->flag2[i] = false;
				quantizer->flag3[i] = false;
				quantizer->flag4[i] = false;
				quantizer->x[i][0] = state[i];
				quantizer->x[i][1] = stateDer[i];
				//printf("x[i][0]: %f\n", x[i][0]);
				//printf("x[i][1]: %f\n", x[i][1]);
				quantizer->q[i][0] = state[i];
				quantizer->x[i][2] = 0;
				quantizer->u0[i] = quantizer->x[i][1] - quantizer->q[i][0]*quantizer->a[i];
				quantizer->u1[i] = 2*quantizer->x[i][2] - quantizer->q[i][1]*quantizer->a[i];
				//printf("u0[i]: %f\n", u0[i]);
				//printf("u1[i]: %f\n", u1[i]);

				quantizer->diffxq[i][1] = quantizer->q[i][1] - quantizer->x[i][1];
				quantizer->diffxq[i][2] = -quantizer->x[i][2];
				quantizer->lqu[i] = quantizer->x[i][0] * tolerance;
				if(quantizer->lqu[i] < absTolerance)
					quantizer->lqu[i] = absTolerance;

				quantizer->diffxq[i][0] = quantizer->q[i][0] - quantizer->dq[i] + quantizer->lqu[i] - quantizer->x[i][0];
				mpr1 = minRootPos(quantizer->diffxq,i,2);
				//printf("mpr1: %f\n", mpr1);
				quantizer->diffxq[i][0] = quantizer->q[i][0] - quantizer->dq[i] - quantizer->lqu[i] - quantizer->x[i][0];
				mpr2 = minRootPos(quantizer->diffxq,i,2);
				//printf("mpr2: %f\n", mpr2);

				if(mpr1 < mpr2)
					quantizer->mpr[i] = mpr1;
				else
					quantizer->mpr[i] = mpr2;
				if(quantizer->mpr[i] > ft)
					quantizer->mpr[i] = ft;
				}
			dt = quantizer->mpr[0];
			dt_index = 0;
			for(i=0;i<STATES;i++){
				if(dt > quantizer->mpr[i]){
					dt = quantizer->mpr[i];
					dt_index = i;
				}
			}
			t=t+dt;
		}
		else{
			printf("else\n");

			if(LIQSS_write(quantizer, data, 3, t, STATES, dt_index) == 1){
				return 1;
			}

			LIQSS2_recomputeNextTimes(quantizer, t, ft, dt_index);

			for(j=0;j<quantizer->nSD[dt_index];j++){
				i= quantizer->SD[dt_index][j];
			}

			t = 1.0 / 0.0;
			for(j=0;j<quantizer->nSD[dt_index];j++){
				i= quantizer->SD[dt_index][j];
				if(quantizer->nTime[i] < t){
					t = quantizer->nTime[i];
					dt_index = i;
				}
			}
			solverInfo->currentTime = t;
		}
		TRACE_POP
	}
    printf("End of simulation, yo!\n");
	sim_result.emit(&sim_result, data, threadData);
	LIQSS2_freeQuantizer (quantizer, STATES);
	TRACE_POP /* pop loop */
	return 0;
}

// For bball on Flat surface
static modelica_real ddx(modelica_real** x, const uinteger index, modelica_real* state){
	if(index == 0){
		return 0.5*(-0.1*x[0][1]);
	}
	else if(index == 1){
		return 0.5*(-0.1*x[1][1] - state[6]*(x[2][1]*1e6 + 30*x[1][1]));
	}
	else if(index == 2){
		return 0.5*(x[1][1]);
	}
	else
		exit(2);
}
