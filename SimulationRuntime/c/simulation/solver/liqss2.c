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
// #include "simulation/solver/liqss2.h"
#include <unistd.h>

// Debugging flag

const bool LIQSS2_DEBUG = false;


//
// void * checkedMalloc (unsigned long long len)
// {
//   void *p = malloc (len);
//   if (!p)
//     {
//       fprintf (stderr, "\nRan out of memory!\n");
//       exit (1);
//     }
//   return ((p));
// }


// LIQSS2_quantizerState LIQSS2_QuantizerState ();









/*
 *  \param [ref] [data]
 *  \param [ref] [solverInfo]
 *
 *  This function performs the simulation controlled by solverInfo.
 */
int prefixedName_LIQSS2Simulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo){
	TRACE_PUSH
	printf("Inside the LIQSS2 solver\n");

	// We need to calculate this Jacobian to ensure that we can use nSD and SD
	data->callback->initialAnalyticJacobianA(data, threadData);

    // Initialization of all the necessary global variables.
    const uinteger STATES = data->modelData->nStates; //Number of state variables

	modelica_real ft = data->simulationInfo->stopTime;

	int iterations = 0;
	modelica_real t = 0;
    const uinteger order = 2;
	const SPARSE_PATTERN* pattern = &(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A].sparsePattern);
	modelica_real dt = 0;
	int dt_index = 0;
	modelica_real** x = NULL;
	const modelica_real tolerance = 0.000001;
	const modelica_real absTolerance = 0.000000001;

	LIQSS2_quantizerState sss;
	printf("sss->order: %d\n", sss->order);
	sss = LIQSS2_QuantizerState ();
	// int u01 = sss->order;
	printf("sss->order: %d\n", sss->order);
	LIQSS2_init(sss, STATES);
	// printf("ss->order: %d\n", ss->order);
	//modelica_real *u01 = ss->order;
	//
	// printf("ss->order: %d\n", u01);
	// // ss->order =5;
	// u01 = 5;
	// printf("ss->order: %d\n", u01);
	printf("sss->order: %d\n", sss->order);
	LIQSS2_freeQuantizer(sss, STATES);
	//state LIQSS2_QuantizerState (STATES);

	modelica_real* tx = NULL;
	modelica_real** q = NULL;
	modelica_real* a = NULL;
	modelica_real* mpr =NULL;
	modelica_real* u0 =NULL;
	modelica_real* u1 =NULL;
	modelica_real** diffxq = NULL;
	modelica_real* lqu =NULL;
	modelica_real* qAux =NULL;
	modelica_real* oldDx =NULL;
	modelica_real* lquOld =NULL;
	modelica_real* dx =NULL;
	modelica_boolean* flag2=NULL;
	modelica_boolean* flag3=NULL;
	modelica_boolean* flag4=NULL;
	modelica_real* dq =NULL;
	modelica_real* lt =NULL;
	modelica_real* tq =NULL;
	modelica_real* ltq =NULL;
	modelica_real* nTime =NULL;
	modelica_integer* nSD =NULL;
	modelica_integer** SD =NULL;
	modelica_boolean fail = 0;  /* When the solver fails */
	modelica_real elapsed;

    // Allocating memory
    uinteger i=0, j=0;
    x = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
	q = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
	diffxq = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
	SD = (modelica_integer**) malloc(sizeof(modelica_integer*)*STATES);
    for(i=0; i<STATES; i++){
        x[i] = (modelica_real*) malloc(sizeof(modelica_real)*(order+1));
		q[i] = (modelica_real*) malloc(sizeof(modelica_real)*(order));
		diffxq[i] = (modelica_real*) malloc(sizeof(modelica_real)*(order+1));
		SD[i] = (modelica_integer*) malloc(sizeof(modelica_integer)*(STATES));
	}

	tx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	//fail = (tx == NULL) ? 1 : ( 0 | fail);
	a = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	mpr = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	u0 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	u1 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	lqu = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	qAux = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	oldDx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	lquOld = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	dx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	flag2 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
	flag3 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
	flag4 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
	dq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	lt = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	tq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	nSD = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
	nTime = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	ltq = (modelica_real*)malloc(STATES * sizeof(modelica_real));

	if (fail)
		printf("BIG ISSUE\n");

	// A summary of which variables are in which derivatives.
	uinteger start =0;
	uinteger jIndex = 0;

	for(i=0;i<STATES;i++){
		start=0;
		if(i > 0)
			start = pattern->leadindex[i-1];
		for(uinteger j = start; j < pattern->leadindex[i]; j++){
			printf("%d %d   %d\n",i,j,pattern->index[j]);
		}
		printf("\n");
	}
	printf("\n");

	jIndex=0;
	start=0;

	for(i=0; i<STATES; i++){
		if(i>0)
			start = pattern->leadindex[i-1];
		printf("We are now considering start=%d\n",start);
		nSD[i] = pattern->leadindex[i] - start;
		printf("pattern->leadindex[i]: %d\n", pattern->leadindex[i]);
		printf("start: %d\n", start);
		printf("nSD[i]: %d\n", nSD[i]);
		jIndex = 0;
		for(uinteger j = start; j < pattern->leadindex[i]; j++){
			printf("\t %d %d %d\n",j, jIndex,pattern->index[j]);
			printf("\t %d\n", pattern->index[j]);
			SD[i][jIndex] = pattern->index[j];
			printf("\t %d\n", pattern->index[j]);
			printf("\t %d\n", SD[i][jIndex]);
			printf("\t %d %d %d %d\n",i,j,jIndex,SD[i][jIndex]);
			jIndex++;
		}
	}

	// Initial values and the initial derivative values
	SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
	modelica_real* state = sData->realVars;
    modelica_real* stateDer = sData->realVars + data->modelData->nStates;

	modelica_real QEvent = false;

	//Writing names of state variables

	//Create the .csv files
	char buf[100];
	for(i=0;i<STATES;i++){
		sprintf(buf, "%s.txt", "/home/anton/Desktop/Results/StateVars");
		FILE *fs;
		if(i==0)
			fs = fopen(buf, "w");
		else
			fs = fopen(buf, "a");
		if(fs == NULL){
			printf("Couldn't open file\n");
			return 1;
		}
		fprintf(fs,"%s\n", data->modelData->realVarsData[i].info.name);
		fclose(fs);
	}

	//printing all the stateVar names
	printf("State var names:\n");
	for(i=0;i< data->modelData->nVariablesReal;i++){
		printf("%s\n",data->modelData->realVarsData[i].info.name);
	}
	printf("The others:\n");

	for(i = 0; i < data->modelData->nVariablesInteger; i++){
		printf("tttt\n");
		printf("%s\n", data->modelData->integerVarsData[i].info.name);
	}
	modelica_real constant;

	for(i = 0; i < data->modelData->nVariablesBoolean; i++){
		printf("pppp\n");
		printf("%s\n", data->modelData->booleanVarsData[i].info.name);
	}


	for(i = 0; i < data->modelData->nVariablesInteger; i++){
		printf("wwww\n");
		printf("%s\n", data->modelData->integerVarsData[i].info.name);
	}

	for(i = 0; i < data->modelData->nVariablesString; i++){
		printf("rrrr\n");
		printf("%s\n", data->modelData->stringVarsData[i].info.name);
	}


	for(i = 0; i < data->modelData->nParametersInteger; i++){
		printf("1111\n");
		printf("%s\n", data->modelData->integerParameterData[i].info.name);
	}

	for(i = 0; i < data->modelData->nParametersReal; i++){
		printf("22222\n");
		printf("%s\n", data->modelData->realParameterData[i].info.name);
	}

	printf("vals\n");
	for(i=0;i< data->modelData->nVariablesReal;i++){
 	   printf("state: %f\n",state[i]);
    }

	//Initialising the algorithm
	modelica_real mpr1 = 0;
	modelica_real mpr2 = 0;
	for(i=0;i<STATES;i++){
		tx[i] = 0;
		tq[i] = 0;
		ltq[i] = 0;
		lt[i] = 0;
		flag2[i] = false;
		flag3[i] = false;
		flag4[i] = false;
		x[i][0] = state[i];
		x[i][1] = stateDer[i];


		// printf("x[i][0]: %f\n", x[i][0]);
		// printf("x[i][1]: %f\n", x[i][1]);
		q[i][0] = state[i];
		x[i][2] = 0;
		u0[i] = x[i][1] - q[i][0]*a[i];
		u1[i] = 2*x[i][2] - q[i][1]*a[i];
		// printf("u0[i]: %f\n", u0[i]);
		// printf("u1[i]: %f\n", u1[i]);

		diffxq[i][1] = q[i][1] - x[i][1];
		diffxq[i][2] = -x[i][2];
		lqu[i] = x[i][0] * tolerance;
		if(lqu[i] < absTolerance)
			lqu[i] = absTolerance;

		diffxq[i][0] = q[i][0] - dq[i] + lqu[i] - x[i][0];
		mpr1 = minRootPos(diffxq,i,2);
		//printf("mpr1: %f\n", mpr1);
		diffxq[i][0] = q[i][0] - dq[i] - lqu[i] - x[i][0];
		mpr2 = minRootPos(diffxq,i,2);
		//printf("mpr2: %f\n", mpr2);

		if(mpr1 < mpr2)
			mpr[i] = mpr1;
		else
			mpr[i] = mpr2;
		if(mpr[i] > ft)
			mpr[i] = ft;
	}

	//Create the .csv files
	for(i=0;i<STATES;i++){
		sprintf(buf, "%s%s.csv", "/home/anton/Desktop/Results/", data->modelData->realVarsData[i].info.name);
		FILE *fs = fopen(buf, "w");
		if(fs == NULL){
			printf("Couldn't open file\n");
			return 1;
		}
		fprintf(fs,"time,q,dq,ddq,x,dx,ddx\n");
		fclose(fs);
	}

	for(i=0;i<STATES;i++){
		sprintf(buf, "%s%s.csv", "/home/anton/Desktop/Results/", data->modelData->realVarsData[i].info.name);
		FILE *fs = fopen(buf, "a");
		if(fs == NULL){
			printf("Couldn't open file\n");
			return 1;
		}
		fprintf(fs,"%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",t,q[i][0],q[i][1],q[i][2],x[i][0],x[i][1],x[i][2]);
		fclose(fs);
	}
	printf("\nThe Initialization stage\n");
	dt = mpr[0];
	dt_index = 0;
	for(i=0;i<STATES;i++){
		printf("nTime[i]:  %d  %.16f\n", i, mpr[i]);
		if(dt > mpr[i]){
			dt = mpr[i];
			dt_index = i;
		}
	}
	t=dt;
	solverInfo->currentTime = t;
	// printf("dt_index: %d\n", dt_index);
	// printf("dt: %f\n", dt);
	// printf("tx[dt_index]: %f\n", tx[dt_index]);

	//Now we do the integration iteratively
	printf("Hier is die simulasie\n");
	while (t < ft){
		iterations += 1;
		printf("\nBeginning of iteration\n");
		printf("t: %.16f\n", t);
		printf("dt_index: %d\n", dt_index);
		// printf("t - tx[dt_index]: %f\n", t - tx[dt_index]);
		//printf("x[dt_index][0]: %f\n", x[dt_index][0]);
		//printf("x[dt_index][1]: %f\n", x[dt_index][1]);
		//printf("x[dt_index][2]: %f\n", x[dt_index][2]);
		elapsed = t - tx[dt_index];
		printf("Elapsed: %.16f\n", elapsed);
		x[dt_index][0] = calcState(x, dt_index, elapsed);
		// printf("x[dt_index][0]: %.16f\n", x[dt_index][0]);

		x[dt_index][1] = calcDerivative(x, dt_index, elapsed);
		// printf("x[dt_index][1]: %.16f\n", x[dt_index][1]);

		tx[dt_index] = t;
		lqu[dt_index] = fabs(x[dt_index][0]) * tolerance;
		if(lqu[dt_index] < absTolerance){
			lqu[dt_index] = absTolerance;
		}
		// printf("lqu[dt_index]: %.12f\n",lqu[dt_index]);

		// printf("Voor QA_updateQuantizedState\n");
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

		// Hier begin QA_updateQuantizedState
		// printf("QA_updateQuantizedState\n");

		flag3[dt_index] = false;
		elapsed = t - tq[dt_index];
		qAux[dt_index] = q[dt_index][0] + elapsed * q[dt_index][1];
		oldDx[dt_index] = x[dt_index][1];
		elapsed = t - lt[dt_index];
		ltq[dt_index] = t;
		u0[dt_index] = u0[dt_index] + elapsed * u1[dt_index];

		// printf("flag2[dt_index] %d\n", flag2[dt_index]);

		if(flag2[dt_index]){
			lqu[dt_index] = lquOld[dt_index];
			flag2[dt_index] = false;
			q[dt_index][0] = qAux[dt_index];
		}
		else
			q[dt_index][0] = x[dt_index][0];

		// printf("!!q[dt_index][0]: %f\n",q[dt_index][0]);
		if(a[dt_index] < -1e-30){
			if(x[dt_index][2] < 0){
				dx[dt_index] = a[dt_index] * a[dt_index] * (q[dt_index][0] + lqu[dt_index]) + a[dt_index] * u0[dt_index] + u1[dt_index];
				if(dx[dt_index] <= 0)
					dq[dt_index] = lqu[dt_index];
				else{
					dq[dt_index] = (-u1[dt_index] / a[dt_index] / a[dt_index]) - (u0[dt_index] / a[dt_index]) - q[dt_index][0];
					flag3[dt_index] = true;
					if(fabs(dq[dt_index]) > lqu[dt_index])
						dq[dt_index] = lqu[dt_index];
				}
			}
			else{
				dx[dt_index] = a[dt_index] * a[dt_index] * (q[dt_index][0] - lqu[dt_index]) + a[dt_index] * u0[dt_index] + u1[dt_index];
				if(dx[dt_index] >= 0)
					dq[dt_index] = -lqu[dt_index];
				else{
					dq[dt_index] = (-u1[dt_index] / a[dt_index] / a[dt_index]) - (u0[dt_index] / a[dt_index]) - q[dt_index][0];
					flag3[dt_index] = true;
					if(fabs(dq[dt_index]) > lqu[dt_index]){
						dq[dt_index] = -lqu[dt_index];
					}
				}
			}
			if(q[dt_index][1] * x[dt_index][1] < 0 && !flag2[dt_index] && !flag3[dt_index] && !flag4[dt_index]){
				if(q[dt_index][1] < 0)
					dq[dt_index] = qAux[dt_index] - q[dt_index][0] - fabs(lquOld[dt_index]) * 0.1;
				else
					dq[dt_index] = qAux[dt_index] - q[dt_index][0] + fabs(lquOld[dt_index]) * 0.1;
				flag4[dt_index] = true;
			}
			else if (flag4[dt_index]){
				flag4[dt_index] = false;
				if(fabs((-u1[dt_index] / a[dt_index] / a[dt_index]) - (u0[dt_index] / a[dt_index]) - q[dt_index][0]) < 3* lqu[dt_index]){
					dq[dt_index] = (-u1[dt_index] / a[dt_index] / a[dt_index]) - (u0[dt_index] / a[dt_index]) - q[dt_index][0];
					flag3[dt_index] = true;
				}
			}
		}
		else{
			flag4[dt_index] = false;
			if(x[dt_index][2] < 0)
				dq[dt_index] = -lqu[dt_index];
			else
				dq[dt_index] = lqu[dt_index];
		}


		// printf("22dq[dt_index]: %f\n",dq[dt_index]);
		if(fabs(dq[dt_index]) > 2* lqu[dt_index]){
			if(dq[dt_index] > 0)
				dq[dt_index] = lqu[dt_index];
			else
				dq[dt_index] = -lqu[dt_index];
		}

		q[dt_index][0] = q[dt_index][0] + dq[dt_index];
		// printf("!!11dq[dt_index]: %f\n",dq[dt_index]);
		// printf("!!11q[dt_index][0]: %f\n",q[dt_index][0]);
		if(flag3[dt_index])
			q[dt_index][1] = a[dt_index] * q[dt_index][0] + u0[dt_index];
		else
			q[dt_index][1] = x[dt_index][1];

		// Hier einding QA_updateQuantizedState -- liqss2.c
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

		tq[dt_index] = t;

		// printf("QA_nextTime\n");
		// Hier begin QA_nextTime -- liqss2.c
		if(x[dt_index][2] == 0)
			nTime[dt_index] = 1.0 / 0.0;
		else
			nTime[dt_index]  = t + sqrt(fabs(lqu[dt_index] / x[dt_index][2]));
		// Hier eindig QA_nextTime -- liqss2.c
		// printf("After QA_nextTime\n");
		for(i=0;i<nSD[dt_index];i++){

			j= SD[dt_index][i];
			printf("j: %d\n", j);
			elapsed = t - tx[j];
			if(elapsed > 0){
				x[j][0] = calcState(x,j,elapsed);
				// printf("x[j][0]: %.16f\n", x[j][0]);
				tx[j] = t;
			}

			// Hierdie is FRW_recomputeDerivatives
			elapsed = t - tq[j];
			if(elapsed > 0){
				q[j][0] = calcQState(q,j,elapsed);
				//printf("q[j][0]: %d   %f\n",j,q[j][0]);
				tq[j] = t;
			}
			// Hier eindig FRW_recomputeDerivatives
		}

		// for(i = 0;i<STATES;i++){
		// 	for(j =0;j<order+1;j++){
		// 		printf("q[i][j]: %d %d   %.16f\n",i,j,q[i][j]);
		// 	}
		// }


		SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
	    modelica_real* state = sData->realVars;
	    modelica_real* stateDer = sData->realVars + data->modelData->nStates;
		uinteger i = 0;
		for(i = 0; i< data->modelData->nStates; i++){
			state[i] = q[i][0];
			// printf("state[i]: %d  %f\n",i,state[i]);
		}
		sData->timeValue = solverInfo->currentTime;
		QEvent = LIQSS2_calculateState(data, threadData);
		simulationUpdate(data, threadData, solverInfo);
		sim_result.emit(&sim_result, data, threadData);
		for(i = 0; i< data->modelData->nStates; i++){
			x[i][1] = stateDer[i];
			// printf("x[i][1]: %d  %f\n",i,x[i][1]);
			// printf("sdsdsdds\n");
		}
		for(i = 0; i< STATES; i++){
			x[i][2] = ddx(q,i,state);
			// printf("x[i][2]: %d  %.16f\n",i,x[i][2]);
		}
		// for(i = 0;i<STATES;i++){
		// 	for(j =0;j<order+1;j++){
		// 		printf("x[i][j]: %d %d   %.16f\n",i,j,x[i][j]);
		// 	}
		// }

		if(QEvent){
			for(i=0;i<STATES;i++){
				tx[i] = t;
				tq[i] = t;
				ltq[i] = t;
				lt[i] = t;
				flag2[i] = false;
				flag3[i] = false;
				flag4[i] = false;
				x[i][0] = state[i];
				x[i][1] = stateDer[i];
				//printf("x[i][0]: %f\n", x[i][0]);
				//printf("x[i][1]: %f\n", x[i][1]);
				q[i][0] = state[i];
				x[i][2] = 0;
				u0[i] = x[i][1] - q[i][0]*a[i];
				u1[i] = 2*x[i][2] - q[i][1]*a[i];
				//printf("u0[i]: %f\n", u0[i]);
				//printf("u1[i]: %f\n", u1[i]);

				diffxq[i][1] = q[i][1] - x[i][1];
				diffxq[i][2] = -x[i][2];
				lqu[i] = x[i][0] * tolerance;
				if(lqu[i] < absTolerance)
					lqu[i] = absTolerance;

				diffxq[i][0] = q[i][0] - dq[i] + lqu[i] - x[i][0];
				mpr1 = minRootPos(diffxq,i,2);
				//printf("mpr1: %f\n", mpr1);
				diffxq[i][0] = q[i][0] - dq[i] - lqu[i] - x[i][0];
				mpr2 = minRootPos(diffxq,i,2);
				//printf("mpr2: %f\n", mpr2);

				if(mpr1 < mpr2)
					mpr[i] = mpr1;
				else
					mpr[i] = mpr2;
				if(mpr[i] > ft)
					mpr[i] = ft;
				}
			dt = mpr[0];
			dt_index = 0;
			for(i=0;i<STATES;i++){
				if(dt > mpr[i]){
					dt = mpr[i];
					dt_index = i;
				}
			}
			t=t+dt;
		}
		else{

			// Begin of QA_recomputeNextTimes
			printf("else\n");

			sprintf(buf, "%s%s.csv", "/home/anton/Desktop/Results/", data->modelData->realVarsData[dt_index].info.name);
			FILE *fs = fopen(buf, "a");
			if(fs == NULL){
				printf("Couldn't open file\n");
				return 1;
			}
			fprintf(fs,"%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",t,q[dt_index][0],q[dt_index][1],q[dt_index][2],x[dt_index][0],x[dt_index][1],x[dt_index][2]);
			fclose(fs);

			modelica_real diffQ = 0, timeaux = 0;
			for(j=0;j<nSD[dt_index];j++){
				i= SD[dt_index][j];
				if(ltq[i] == t){
					diffQ = q[i][0] - qAux[i];
					if(fabs(diffQ) > lqu[i] * 1e-6){
						a[i] = (x[i][1] - oldDx[i]) / diffQ;
						if(a[i] > 0)
							a[i] = 0;
					}
				}
				else
					flag3[i] = false;

				u0[i] = x[i][1] - q[i][0] * a[i];
				u1[i] = 2 * x[i][2] - q[i][1] * a[i];
				lt[i] = t;
				if(flag4[i])
					nTime[i] = t;
				else{
					diffxq[i][1] = q[i][1] - x[i][1];
					diffxq[i][2] = -x[i][2];
					diffxq[i][0] = q[i][0] - dq[i] + lqu[i] - x[i][0];
					nTime[i] = t + minRootPos(diffxq,i,2);
					diffxq[i][0] = q[i][0] - dq[i] - lqu[i] - x[i][0];
					timeaux = t + minRootPos(diffxq,i,2);
					if(timeaux < nTime[i])
						nTime[i] = timeaux;
					if(a[i] != 0 && fabs(x[i][2]) > 1e-10 && !flag3[i] && !flag2[i]){
						diffxq[i][0] = a[i] * a[i] * q[i][0] + a[i] * u0[i] + u1[i];
						diffxq[i][1] = a[i] * a[i] * q[i][1] + a[i] * u1[i];
						timeaux = t + minRootPos(diffxq, i, 1);
						if(timeaux < nTime[i]){
							flag2[i] = true;
							nTime[i] = timeaux;
							lquOld[i] = lqu[i];
						}
					}
					else
						flag2[i] = false;

					if(nTime[i] > ft)
						nTime[i] = ft;

					// double err1 = q[i] - x[i] + diffxq[i][1] * (nTime[i] - t) / 2 + diffxq[i][2] * pow((nTime[i] - t) / 2, 2);
					// //double err1 = 0;
					// if(fabs(err1) > 3 * fabs(lqu[i]))
					// 	nTime[i] = t + ft * 1e-3;
				}
			}

			// End of QA_recomputeNextTimes


			for(j=0;j<nSD[dt_index];j++){
				i= SD[dt_index][j];
				printf("nTime[i]: %d  %.16f\n",i,nTime[i]);
			}

			t = 1.0 / 0.0;
			//printf("fdf\n");
			//printf("t: %.12f\n",t);
			for(j=0;j<nSD[dt_index];j++){
				i= SD[dt_index][j];
				if(nTime[i] < t){
					t = nTime[i];
					dt_index = i;
				}
			}
			//printf("t: %.12f\n",t);
			solverInfo->currentTime = t;
			//printf("End of simulation, yo!aaa\n");
			//printf("dt_index: %d\n",dt_index);
		}
		TRACE_POP
	}

    printf("End of simulation, yo!\n");

	sim_result.emit(&sim_result, data, threadData);



	/* free memory*/
    for(i=0; i<STATES; i++){
        free(x[i]);
		free(q[i]);
		free(diffxq[i]);
		//free(SD[i]);
	}
	free(a);
	free(mpr);
	free(u0);
	free(u1);
	free(lqu);
	free(qAux);
	free(oldDx);
	free(lquOld);
	free(dx);
	free(flag2);
	free(flag3);
	free(flag4);
	free(dq);
	free(lt);
	free(tq);
	free(nSD);
	free(nTime);
	free(ltq);
	for(i=0; i<STATES; i++){
		free(SD[i]);
	}
	free(tx);
	/* end - free memory */

	TRACE_POP /* pop loop */

	return 0;
}

// static modelica_boolean calcAllDerivatives(modelica_real** q, modelica_real** x, DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo){
// 	printf("In calcAllDerivatives\n");
// 	SIMULATION_DATA *sData = (SIMULATION_DATA*)data->localData[0];
//     modelica_real* state = sData->realVars;
//     modelica_real* stateDer = sData->realVars + data->modelData->nStates;
// 	uinteger i = 0;
// 	for(i = 0; i< data->modelData->nStates; i++){
// 		state[i] = q[i][0];
// 		printf("state[i]: %d  %f\n",i,state[i]);
// 	}
// 	sData->timeValue = solverInfo->currentTime;
// 	modelica_boolean QEvent = LIQSS2_calculateState(data, threadData);
// 	simulationUpdate(data, threadData, solverInfo);
// 	for(i = 0; i< data->modelData->nStates; i++){
// 		x[i][1] = stateDer[i];
// 		printf("x[i][1]: %d  %f\n",i,x[i][1]);
// 	}
// 	for(i = 0; i< data->modelData->nStates; i++){
// 		x[i][2] = ddx(q,i);
// 	}
// 	*x = x;
// 	return QEvent;
// }

static modelica_boolean LIQSS2_calculateState(DATA* data, threadData_t *threadData){
	externalInputUpdate(data);
	data->callback->input_function(data, threadData);
	data->callback->functionODE(data, threadData);
	data->callback->functionAlgebraics(data, threadData);
	data->callback->output_function(data, threadData);
	data->callback->function_storeDelayed(data, threadData);
	data->callback->functionDAE(data, threadData);

	if(data->simulationInfo->needToIterate){
		//printf("needToIterate is een\n");
		//printf("nTI:  %d\n", data->simulationInfo->needToIterate);
		return true;
	}
	else
		return false;
}

static modelica_real calcState(modelica_real** x, const uinteger index, const modelica_real dt){
	//printf("In calc state\n");
	//printf("dt %f\n", dt);
	return x[index][0] + dt * x[index][1] + dt * dt * x[index][2];
}

static modelica_real calcDerivative(modelica_real** x, const uinteger index, const modelica_real dt){
	return x[index][1] + 2 * dt * x[index][2];
}

static modelica_real calcQState(modelica_real** q, const uinteger index, const modelica_real dt){
	return q[index][0] + dt * q[index][1];
}

// For Lotka-Volterra
// static modelica_real ddx(modelica_real** x, const uinteger index){
// 	if(index == 0){
// 		return 0.5*(1.5*x[0][1] - x[0][0]*x[1][1] - x[0][1]*x[1][0]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(x[0][0]*x[1][1] + x[0][1]*x[1][0] - 7*x[1][1]);
// 	}
// 	else
// 		exit(2);
// }



// static modelica_real ddx(modelica_real** x, const uinteger index, modelica_real* state){
// 	if(index == 0){
// 		return 0.5*(-0.1*x[0][1]);
// 	}
// 	else if(index == 1){
// 		return 0.5*(-0.1*x[1][1] - state[8]*(x[3][1]*1e6 + 30*x[1][1]));
// 	}
// 	else if(index == 2){
// 		return 0.5*(x[0][1]);
// 	}
// 	else if(index == 3){
// 		return 0.5*(x[1][1]);
// 	}
// 	else
// 		exit(2);
// }


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



// static modelica_real ddx(modelica_real** x, const uinteger index, modelica_real* state){
// 	modelica_real mu =1000;
// 	if(index > 0 && index < 500){
// 		return 0.5*(mu*x[index][1]*(0.5-x[index][0])*x[index][0] -(-1.0+x[index][0])*mu*x[index][1]*x[index][0] -(x[index][1]-x[index-1][1])*500+ (-1.0+x[index][0])*mu*x[index][1]*(0.5-x[index][0]));
// 	}
// 	else if(index == 0){
// 		return 0.5*(-500*x[index][1]-x[index][0]*mu*(-1.0+x[index][0])*x[index][1]-(x[index][0]-0.5)*mu*(-1.0+x[index][0])*x[index][1]-x[index][0]*(x[index][0]-0.5)*mu*x[index][1]);
// 	}
// 	else{
// 		printf("Groot probleme hier\n");
// 	}
//
// }


// static modelica_real ddx(modelica_real** x, const uinteger index, modelica_real* state){
// 	return 0.5;
// }






static modelica_real minRootPos(modelica_real** diffxq, const uinteger index, const uinteger order){
	modelica_real mpr = 0;
	if(order == 1){
		if(diffxq[index][1] == 0)
			return 1.0 /0.0;
		else{
			mpr = -diffxq[index][0] / diffxq[index][1];
			if(mpr < 0)
				return 1.0 /0.0;
			else
				return mpr;
		}
	}
	else if(order == 2){
		//printf("diffxq[index][0]: %.12f\n", diffxq[index][0]);
        //printf("diffxq[index][1]: %.12f\n", diffxq[index][1]);
        //printf("diffxq[index][2]: %.12f\n", diffxq[index][2]);
		if(diffxq[index][2] == 0 || (1000*fabs(diffxq[index][2]) < fabs(diffxq[index][1]))){
			if(diffxq[index][1] == 0)
				return 1.0 /0.0;
			else{
				mpr = -diffxq[index][0] / diffxq[index][1];
				if(mpr < 0)
					return 1.0 / 0.0;
				else
					return mpr;
			}
		}
		else{
			modelica_real disc = diffxq[index][1] * diffxq[index][1] - 4 * diffxq[index][2] * diffxq[index][0];
			if(disc < 0)
				return 1.0 /0.0;
			else{
				//printf("disc: %.12f\n", disc);
				modelica_real sd = sqrt(disc);
				//printf("sd: %.12f\n", sd);
				modelica_real r1 = (-diffxq[index][1] + sd) / (2 * diffxq[index][2]);
				if(r1 > 0)
					mpr = r1;
				else
					mpr = 1.0/0.0;
				r1 = (-diffxq[index][1] - sd) / (2 * diffxq[index][2]);
				if(r1 > 0 && r1 < mpr){
					mpr = r1;
				}
				//printf("mpr: %.12f\n",mpr);
				return mpr;
			}
		}
	}
	else
		exit(1);
}
