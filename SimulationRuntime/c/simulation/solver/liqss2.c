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
	uinteger start =0,jIndex = 0, i = 0, j = 0, k = 0;

	const modelica_real tolerance = 1e-6;
	const modelica_real absTolerance = 1e-9;

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
		//quantizer->SD[i][jIndex] = i;
	}

	// for(i =0;i<STATES;i++){
	// 	printf("nSD[i]  %d  %d   \n", i,quantizer->nSD[i]);
	// 	for(j =0; j<quantizer->nSD[i]; j++){
	// 		printf("SD[i][j]  %d  %d   %d\n", i,j,quantizer->SD[i][j]);
	// 	}
	// 	printf("\n");
	// }


	// Initial values and the initial derivative values
	SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
	modelica_real* state = sData->realVars;
    modelica_real* stateDer = sData->realVars + data->modelData->nStates;
	modelica_real QEvent = false;


	if(LIQSS_write(quantizer, data, 1, t, STATES, dt_index) == 1){
		return 1;
	}

	LIQSS2_printModelData(data);

	LIQSS_initialiation(quantizer, state, stateDer, STATES, ft, tolerance, absTolerance);

	if(LIQSS_write(quantizer, data, 2, t, STATES, dt_index) == 1){
		return 1;
	}


	// printf("\nThe Initialization stage\n");
	dt = quantizer->mpr[0];
	dt_index = 0;
	for(i=0;i<STATES;i++){
		// printf("nTime[i]:  %d  %.20f\n", i, quantizer->mpr[i]);
		if(dt > quantizer->mpr[i]){
			dt = quantizer->mpr[i];
			dt_index = i;
		}
	}
	// printf("dt_index: %d\n", dt_index);
	t=dt;
	solverInfo->currentTime = t;

	for(i=0;i<STATES;i++){
		quantizer->q[i][2] = 0;
		quantizer->nTime[i] = ft;
	}


	// Now we do the integration iteratively
	// This is based on the QSS_SEQH_integrate function in qss_seqh_integrator.c

	while (t < ft){
	// ft=0;
	// while(ft<1){
	// 	ft=22;
		iterations += 1;

		printf("---------------------------------\n");
		elapsed = t - quantizer->tx[dt_index];
		// printf("Elapsed: %.16f\n", elapsed);
		// advanceTime in src/common/utils.c
		quantizer->x[dt_index][0] = calcState(quantizer->x, dt_index, elapsed);
		// printf("elapsed: %.20f\n", elapsed);
		printf("t: %.20f\n", t);
		printf("index: %d\n", dt_index);
		printf("iterations: %d\n", iterations);
		quantizer->x[dt_index][1] = calcDerivative(quantizer->x, dt_index, elapsed);




		quantizer->tx[dt_index] = t;
		quantizer->lqu[dt_index] = fabs(quantizer->x[dt_index][0]) * tolerance;
		if(quantizer->lqu[dt_index] < absTolerance){
			quantizer->lqu[dt_index] = absTolerance;
		}

		// printf("WBEFORE\n");
		// for(i = 0;i <STATES; i++){
		// 	for(j = 0;j<3;j++){
		// 		printf("q[i]: %d   %.20f\n", 3*i+j, quantizer->q[i][j]);
		// 	}
		// }
		// for(i = 0;i <STATES; i++){
		// 	for(j = 0;j<3;j++){
		// 		printf("x[i]: %d   %.20f\n", 3*i +j, quantizer->x[i][j]);
		// 	}
		// }
		// printf("\n");


		// QA_updateQuantizedState -- liqss2.c
		LIQSS2_updateQuantizedState(quantizer, t, dt_index);

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
		}

		for(j=0;j<STATES;j++){
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
			// state[i] = i*123;
			// quantizer->q[i][0] = i*123;
		}
		sData->timeValue = solverInfo->currentTime;
		QEvent = LIQSS2_calculateState(data, threadData);
		simulationUpdate(data, threadData, solverInfo);
		sim_result.emit(&sim_result, data, threadData);

		for(j=0;j<quantizer->nSD[dt_index];j++){
			i = quantizer->SD[dt_index][j];
		// for(i=0;i<STATES;i++){
			quantizer->x[i][1] = stateDer[i];
			// quantizer->x[i][1] = dx(data, t,quantizer->q,i,state);
			// printf("der[%d + 1]:  %.20f\n", i*3, quantizer->x[i][1]);
			// printf("stateDer[%d]:  %.20f\n\n", i, stateDer[i]);
			quantizer->x[i][2] = ddx(data,t,quantizer->q,i,state);
			// printf("der[%d + 2]:  %.20f\n", i*3, quantizer->x[i][2]);
		}
		// FRW_recomputeDerivatives ends here -- SYM_recomputeDerivatives qss_frw_imp.c

		// FRW_nextEventTime -- SYM_nextEventTime qss_seqh_integrator.c
		// There is alot of hackery here for now
		// This code only pertains to the bouncing ball at this stage
		// if (dt_index == 2)
		// 	FRW_nextEventTime(quantizer, t, dt_index);

		if(LIQSS_write(quantizer, data, 3, t, STATES, dt_index) == 1){
			return 1;
		}

		LIQSS2_recomputeNextTimes(STATES, quantizer, t, ft, dt_index);

		t = 1.0 / 0.0;
		int dt_index_temp = 0;
		for(i = 0;i <STATES; i++){
			// printf("nTime[j]:  %d  %.20f\n", i,quantizer->nTime[i]);
			if(quantizer->nTime[i] < t){
				t = quantizer->nTime[i];
				dt_index_temp = i;
			}
		}
		solverInfo->currentTime = t;
		dt_index = dt_index_temp;

		// printf("index: %d\n\n", dt_index);

		TRACE_POP
	}

    printf("End of simulation!\n");
	sim_result.emit(&sim_result, data, threadData);
	LIQSS2_freeQuantizer (quantizer, STATES);
	TRACE_POP /* pop loop */
	return 0;
}




// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(-0.1*x[0][1] - state[4]*(x[1][1]*1e6 + 30*x[0][1]));
// 	}
// 	else if(index == 1){
// 		return 0.5*(x[0][1]);
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }
//
//
//
// // For bball on Flat surface
// static modelica_real dx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 1*(-9.8 -0.1*x[0][0] - state[4]*((x[1][0]-10)*1e6 + 30*x[0][0]));
// 	}
// 	else if(index == 1){
// 		return 1*(x[0][0]);
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }



// static modelica_real dx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	modelica_real __C = 1.000000000000000047921736e-04;
// 	modelica_real __L = 1.000000000000000047921736e-04;
// 	modelica_real __R = 10.0;
// 	modelica_real __U = 24.0;
// 	modelica_real __T = 1.000000000000000047921736e-04;
// 	modelica_real __DC = 2.500000000000000000000000e-01;
// 	modelica_real __ROff = 1.000000000000000000000000e+05;
// 	modelica_real __C1 = 1.000000000000000047921736e-04;
// 	modelica_real __L1 = 1.000000000000000047921736e-04;
//
// 	modelica_real alg0 = (__L*x[1][0]+__L*x[0][0])/__L1;
// 	modelica_real alg1 = (1.0/(__L1))*(__L*x[1][1]+__L*x[0][1]);
// 	modelica_real alg3 = (state[11]*(x[0][0]+alg0)-x[3][0])/(state[10]+state[11]);
// 	modelica_real alg4 = (state[11]*(x[0][1]+alg1)-x[3][1])*(1.0/(state[11]+state[10]));
// 	if(index == 0){
// 		return ((-x[2][0] - alg3*state[10])/(__L));
// 	}
// 	else if(index == 1){
// 		return ((__U+x[2][0]-x[3][0])/(__L));
// 	}
// 	else if(index == 2){
// 		return ((x[0][0]-(x[2][0]/(__R)))/__C);
// 	}
// 	else if(index == 3){
// 		return (alg3 - x[0][0])/__C1;
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }
//
//
// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	modelica_real __C = 1.000000000000000047921736e-04;
// 	modelica_real __L = 1.000000000000000047921736e-04;
// 	modelica_real __R = 10.0;
// 	modelica_real __U = 24.0;
// 	modelica_real __T = 1.000000000000000047921736e-04;
// 	modelica_real __DC = 2.500000000000000000000000e-01;
// 	modelica_real __ROff = 1.000000000000000000000000e+05;
// 	modelica_real __C1 = 1.000000000000000047921736e-04;
// 	modelica_real __L1 = 1.000000000000000047921736e-04;
//
// 	modelica_real alg0 = (__L*x[1][0]+__L*x[0][0])/__L1;
// 	modelica_real alg1 = (1.0/(__L1))*(__L*x[1][1]+__L*x[0][1]);
// 	modelica_real alg3 = (state[11]*(x[0][0]+alg0)-x[3][0])/(state[10]+state[11]);
// 	modelica_real alg4 = (state[11]*(x[0][1]+alg1)-x[3][1])*(1.0/(state[11]+state[10]));
//
//
// 	if(index == 0){
// 		return 0.5*((-x[2][1] - alg4*state[10])*(1.0/(__L)));
// 	}
// 	else if(index == 1){
// 		return 0.5*((x[2][1]-x[3][1])*(1.0/(__L)));
// 	}
// 	else if(index == 2){
// 		return ((x[0][1]-(x[2][1]*(1.0/(__R))))*(1.0/(__C)))/2;
// 	}
// 	else if(index == 3){
// 		return ((1.0/(__C1))*(alg4 - x[0][1]))/2;
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }






// For bball on Flat surface
// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(-0.1*x[0][1]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(-0.1*x[1][1] - state[6]*(x[2][1]*1e6 + 30*x[1][1]));
// 	}
// 	else if(index == 2){
// 		return 0.5*(x[1][1]);
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }

//
//
// // For bball on Flat surface
// static modelica_real dx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 1*(-0.1*x[0][0]);
// 	}
// 	else if(index == 1){
// 		return 1*(-9.8 -0.1*x[1][0] - state[4]*((x[2][0]-10)*1e6 + 30*x[1][0]));
// 	}
// 	else if(index == 2){
// 		return 1*(x[1][0]);
// 	}
// 	else{
// 		printf("BIG PROBLEMS\n");
// 		exit(2);
// 	}
// }


// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.00;
// 	//exit(2);
// }


// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(x[1][1]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(-x[0][1]);
// 	}
// 	else if(index == 2){
// 		return 0.5*(1.5*x[2][1] - x[2][0]*x[3][1] - x[2][1]*x[3][0]);
// 	}
// 	else if(index == 3){
// 		return 0.5*(x[2][0]*x[3][1] + x[2][1]*x[3][0] - 2*x[3][1]);
// 	}
// 	else
// 		exit(2);
// }


// static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.3;
// }
//
// static modelica_real dx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.3;
// }


// Van der Pol
// static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(x[1][1]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(data->simulationInfo->realParameter[0]*x[1][1] - 2*data->simulationInfo->realParameter[0]*x[0][0]*x[0][1]*x[1][0] - data->simulationInfo->realParameter[0]*x[0][0]*x[0][0]*x[1][1] - x[0][1]);
// 		//return 2;
// 	}
// 	else
// 		exit(2);
// }



// static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger i, modelica_real* state){
// 	modelica_real mu = 1000;
// 	modelica_real N = 900;
// 	modelica_real alpha  = 0.5;
// 	if(i > 0 && i < 900){
// 		return 0.5*( -x[i][0]*x[i][1]*(x[i][0]-1.0)*mu - (x[i][0]-0.5)*x[i][0]*x[i][1]*mu + (x[i-1][1] - x[i][1])*900 - (x[i][0] - 0.5)*x[i][1]*(x[i][0] - 1.0)*1000);
// 	}
// 	else if(i == 0){
// 		return 0.5*( -x[i][0]*x[i][1]*(x[i][0]-1.0)*mu - (x[i][0]-0.5)*x[i][0]*x[i][1]*mu - (x[i][1])*900 - (x[i][0] - 0.5)*x[i][1]*(x[i][0] - 1.0)*1000);
// 	}
// 	else{
// 		printf("Groot probleme hier\n");
// 	}
// }



// static modelica_real ddx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(1.5*x[0][1] - x[0][0]*x[1][1] - x[0][1]*x[1][0]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(x[0][0]*x[1][1] + x[0][1]*x[1][0] - 2*x[1][1]);
// 	}
// 	else
// 		exit(2);
// }
//
// static modelica_real dx(modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 1*(1.5*x[0][0] - x[0][0]*x[1][0]);
// 	}
// 	else if(index == 1){
// 		return 1*(x[0][0]*x[1][0] - 1*x[1][0]);
// 	}
// 	else
// 		exit(2);
// }






static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
	modelica_real pecPhi = state[18] * sin(3.141592653589793 * ((t - state[19])*(1/0.45)));
	modelica_real pecE = 0.1 * (1 - pecPhi) + 1.1 * pecPhi;
	if(index == 0){
		return 0.5*(x[1][1] - x[2][1]);
	}
	else if(index == 1){
		return 0.5*(pecE * x[3][1] - 0.1 * x[0][1] + 10.0 * x[1][1]) * (1/0.0001);
	}
	else if(index == 2){
		return 0.5*(pecE * x[0][1] - 0.1 * x[3][1] + 10.0 * x[2][0]) * (1/0.0001);
	}
	else if(index == 3){
		return 0.5 * (x[2][1] - x[1][1]);
	}
	else
		exit(2);
}

// static modelica_real dx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.3;
// }



// static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.3;
// }

// static modelica_real dx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.3;
// }


// static modelica_real ddx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(x[13][1] - x[1][1]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(((0.075 * (x[0][1]) - x[1][1] * 0.000089) + (0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[2][1]))  * (1.0 /  0.00005));
// 		// return 0.5*((((0.075 * (x[0][1]) - x[1][1] * 0.000089) + (0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[2][1]))));
// 		//return 0.3;
// 	}
// 	else if(index == 2){
// 		return 0.5*(x[1][1] - x[3][1]);
// 	}
// 	else if(index == 3){
// 		//return 0.5*((((x[1][1] - x[3][1]))  - (1.2870012870012870 * x[16][1] + 0.0334 * x[3][1])));
// 		return 0.5*((((0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[1][1] - x[3][1]))  - (1.2870012870012870 * x[16][1] + 0.0334 * x[3][1])) * (1.0 / 0.000416));
// 		// return 0.5*((((0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[1][1] - x[3][1]))  - (1.2870012870012870 * x[16][1] + 0.0334 * x[3][1])));
// 		//return 0.3;
// 	}
// 	else if(index == 4){
// 		return 0.5*(x[20][1] - x[5][1]);
// 	}
// 	else if(index == 5){
// 		return 0.5*(((0.06 * x[4][1] - 0.0000594 * x[5][1]) - ((0.0243 * (1 - state[70]) + 0.523 * state[70]) * x[6][1] ))  * (1.0 /  0.00005));
// 		// return 0.5*(((0.06 * x[4][1] - 0.0000594 * x[5][1]) - ((0.0243 * (1 - state[70]) + 0.523 * state[70]) * x[6][1] )));
// 		//return 0.3;
// 	}
// 	else if(index == 6){
// 		return 0.5*(x[5][1] - x[7][1]);
// 	}
// 	else if(index == 7){
// 		return 0.5*((((0.0243 * (1 - state[70]) + 0.523 * state[70]) * x[6][1]) - (0.4500450045004500 * x[9][1] + 0.0251 * x[7][1]))  * (1.0 /  0.000206));
// 		// return 0.5*((((0.0243 * (1 - state[70]) + 0.523 * state[70]) * x[6][1]) - (0.4500450045004500 * x[9][1] + 0.0251 * x[7][1])));
// 		//return 0.3;
// 	}
// 	else if(index == 8){
// 		return 0.5*((0.5624296962879640 * x[11][1] - 0.6752194463200539 * x[8][1]) * (1.0 /  0.053) + x[10][1]);
// 	}
// 	else if(index == 9){
// 		return 0.5*(x[7][1] - x[10][1]);
// 	}
// 	else if(index == 10){
// 		return 0.5*((0.4500450045004500 * x[9][1] - 0.0227 * x[10][1] - 0.6752194463200539 * x[8][1]) * (1.0 /  0.00005));
// 		// return 0.5*((0.4500450045004500 * x[9][1] - 0.0227 * x[10][1] - 0.6752194463200539 * x[8][1]));
// 	}
// 	else if(index == 11){
// 		return 0.5*(((0.1500150015001500 * x[14][1] - 0.5624296962879640 * x[11][1]) * (1.0 /  0.0379)) - ((0.5624296962879640 * x[11][1] - 0.6752194463200539 * x[8][1])  * (1.0 /  0.053)));
// 		// return 0.5*(((0.1500150015001500 * x[14][1] - 0.5624296962879640 * x[11][1]) / 0.0379) - ((0.5624296962879640 * x[11][1] - 0.6752194463200539 * x[8][1])));
// 	}
// 	else if(index == 12){
// 		return 0.5*(- x[13][1] - ((0.2 * x[12][1] - 0.1500150015001500 * x[14][1])  * (1.0 /  0.0252)));
// 		// return 0.5*(- x[13][1] - ((0.2 * x[12][1] - 0.1500150015001500 * x[14][1])));
// 	}
// 	else if(index == 13){
// 		return 0.5*((0.2 * (x[12][1]) - 0.0126 * x[13][1] - 0.075 * (x[0][1])) * (1.0 /0.00005));
// 		// return 0.5*((0.1500150015001500 * x[14][1] - 0.075 * x[0][1]));
// 	}
// 	else if(index == 14){
// 		return 0.5*(((0.2 * x[12][1] - 0.1500150015001500 * x[14][1])  * (1.0 /  0.0252)) - ((0.1500150015001500 * x[14][1] - 0.5624296962879640 * x[11][1])  * (1.0 /  0.0379)));
// 	}
// 	else if(index == 15){
// 		return 0.5*(((0.5524861878453039 * x[18][1] - 0.6097560975609756 * x[15][1])  * (1.0 / 0.178)) + x[17][1]);
// 	}
// 	else if(index == 16){
// 		return 0.5*(x[3][1] - x[17][1]);
// 	}
// 	else if(index == 17){
// 		return 0.5*((1.2870012870012870 * x[16][1] - 0.0824 * x[17][1] - 0.6097560975609756 * x[15][1]) * (1.0 / 0.00005));
// 		// return 0.5*((1.2870012870012870 * x[16][1] - 0.0824 * x[17][1] - 0.6097560975609756 * x[15][1]));
// 	}
// 	else if(index == 18){
// 		return 0.5*(((0.0755287009063444 * x[21][1] - 0.5524861878453039 * x[18][1]) * (1.0 / 0.667)) - ((0.5524861878453039 * x[18][1] - 0.6097560975609756 * x[15][1]) * (1.0 / 0.178)));
// 	}
// 	else if(index == 19){
// 		return 0.5*(- x[20][1] - ((0.0135354629128316 * x[19][1] - 0.0755287009063444 * x[21][1]) *(1.0/0.0233)));
// 		// return 0.5*(- x[20][1] - ((0.0135354629128316 * x[19][1] - 0.0755287009063444 * x[21][1])));
// 	}
// 	else if(index == 20){
// 		return 0.5*((0.0135354629128316 * x[19][1] - 0.0267 * x[20][1] - 0.06 * x[4][1])*(1.0/0.00005));
// 		// return 0.5*((0.0135354629128316 * x[19][1] - 0.0267 * x[20][1] - 0.06 * x[4][1]));
// 	}
// 	else if(index == 21){
// 		return 0.5*((0.0135354629128316 * x[19][1] - 0.0755287009063444 * x[21][1])*(1.0/0.0233) - (0.0755287009063444 * x[21][1] - 0.5524861878453039 * x[18][1])*(1.0/0.667));
// 	}
// 	else
// 		exit(2);
// }
//
//
//
//
//
//
// static modelica_real dx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return (x[13][0] - x[1][0]);
// 	}
// 	else if(index == 1){
// 		// printf("!!!state[55]:  %.20f\n",state[55]);
// 		// printf("!!!$Pheart$PleftHeart$Patrium$Pjunction$Pinductor$Pflow_in$PP:  %.20f\n", (0.075 * (x[0][0] - 30.0) - x[1][0] * 0.000089));
// 		// printf("!!!$Pheart$PleftHeart$PatrioventricularValve$Pflow_in$PP:  %.20f\n", (0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[2][0] - 10.0));
// 		// printf("!!!$Pheart$PleftHeart$Patrium$Pjunction$Pinductor$PL:  0.00005\n");
// 		return ((((0.075 * (x[0][0] - 30.0) - x[1][0] * 0.000089) - (0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[2][0] - 10.0))) / 0.00005);
// 		//return 0.3;
// 	}
// 	else if(index == 2){
// 		return (x[1][0] - x[3][0]);
// 	}
// 	else if(index == 3){
// 		return ((((0.049 * (1 - state[55]) + 2.49 * state[55]) * (x[1][0] - x[3][0] - 10.0))  - (1.2870012870012870 * (x[16][0] - 205.0) + 0.0334 * x[3][0]))/0.000416);
// 		//return 0.3;
// 	}
// 	else if(index == 4){
// 		return (x[20][0] - x[5][0]);
// 	}
// 	else if(index == 5){
// 		return (((0.06 * (x[4][0] - 30.0) - 0.0000594 * x[5][0]) - ((0.0243 * (1 - state[70]) + 0.523 * state[70]) * (x[6][0] - 10.0 ))) / 0.00005);
// 		//return 0.3;
// 	}
// 	else if(index == 6){
// 		return (x[5][0] - x[7][0]);
// 	}
// 	else if(index == 7){
// 		return ((((0.0243 * (1 - state[70]) + 0.523 * state[70]) * (x[6][0] - 10.0)) - (0.4500450045004500 * (x[9][0] - 50.0) + 0.0251 * x[7][0])) / 0.000206);
// 		//return 0.3;
// 	}
// 	else if(index == 8){
// 		return (((0.5624296962879640 * (x[11][0] - 53.0) - 0.6752194463200539 * (x[8][0] - 30.0)) / 0.053) + x[10][0]);
// 	}
// 	else if(index == 9){
// 		return (x[7][0] - x[10][0]);
// 	}
// 	else if(index == 10){
// 		return ((0.4500450045004500 * (x[9][0] - 50.0) - 0.0227 * x[10][0] - 0.6752194463200539 * (x[8][0] - 30.0)) / 0.00005);
// 	}
// 	else if(index == 11){
// 		return (((0.1500150015001500 * (x[14][0] - 75.0) - 0.5624296962879640 * (x[11][0] - 53.0)) / 0.0379) - ((0.5624296962879640 * (x[11][0] - 53.0) - 0.6752194463200539 * (x[8][0] - 30.0)) / 0.053));
// 		}
// 	else if(index == 12){
// 		return (- x[13][0] - ((0.2 * (x[12][0] - 75.0) - 0.1500150015001500 * (x[14][0] - 75.0)) / 0.0252));
// 	}
// 	else if(index == 13){
// 		return (0.2 * (x[12][0] - 75.0) - 0.0126 * x[13][0] - 0.075 * (x[0][0] - 30.0)) / 0.00005;
// 	}
// 	else if(index == 14){
// 		return (((0.2 * (x[12][0] - 75.0) - 0.1500150015001500 * (x[14][0]- 75.0)) / 0.0252) - ((0.1500150015001500 * (x[14][0] - 75.0) - 0.5624296962879640 * (x[11][0] - 53.0)) / 0.0379));
// 	}
// 	else if(index == 15){
// 		return (((0.5524861878453039 * (x[18][0] - 401.0) - 0.6097560975609756 * (x[15][0] - 370.0)) / 0.178) + x[17][0]);
// 	}
// 	else if(index == 16){
// 		return (x[3][0] - x[17][0]);
// 	}
// 	else if(index == 17){
// 		return ((1.2870012870012870 * (x[16][0] - 205.0) - 0.0824 * x[17][0] - 0.6097560975609756 * (x[15][0] - 370.0)) / 0.00005);
// 	}
// 	else if(index == 18){
// 		return (((0.0755287009063444 * (x[21][0] - 596.0) - 0.5524861878453039 * (x[18][0] - 401.0)) / 0.667) - ((0.5524861878453039 * (x[18][0] - 401.0) - 0.6097560975609756 * (x[15][0] - 370.0)) / 0.178));
// 	}
// 	else if(index == 19){
// 		return (- x[20][0] - ((0.0135354629128316 * (x[19][0] - 1938.0) - 0.0755287009063444 * (x[21][0] - 596.0)) / 0.0233));
// 	}
// 	else if(index == 20){
// 		return ((0.0135354629128316 * (x[19][0] - 1938.0) - 0.0267 * x[20][0] - 0.06 * (x[4][0] - 30.0)) / 0.00005);
// 	}
// 	else if(index == 21){
// 		return (((0.0135354629128316 * (x[19][0] - 1938.0) - 0.0755287009063444 * (x[21][0] - 596.0)) / 0.0233) - ((0.0755287009063444 * (x[21][0] - 596.0) - 0.5524861878453039 * (x[18][0] - 401.0)) / 0.667));
// 	}
// 	else
// 		exit(2);
// }
