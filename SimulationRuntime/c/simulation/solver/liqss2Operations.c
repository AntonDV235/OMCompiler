#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solver_main.h"
#include "simulation/simulation_runtime.h"
#include "simulation/results/simulation_result.h"
#include "openmodelica_func.h"
#include "util/omc_error.h"
#include "simulation/options.h"


#include <unistd.h>

void * checkedMalloc (unsigned long long len)
{
  void *p = malloc (len);
  if (!p)
    {
      fprintf (stderr, "\nRan out of memory!\n");
      exit (1);
    }
  return ((p));
}


LIQSS2_quantizerState LIQSS2_QuantizerState (){
  LIQSS2_quantizerState p = checkedMalloc (sizeof(*p));
  p->order = 2;
  p->events = 1;
  // printf("p->order: %d\n", p->order);
  p->x = NULL;
  p->diffxq = NULL;
  p->q = NULL;
  p->SD =NULL;
  p->tx = NULL;
  p->a = NULL;
  p->mpr =NULL;
  p->u0 =NULL;
  p->u1 =NULL;
  p->lqu =NULL;
  p->qAux =NULL;
  p->oldDx =NULL;
  p->lquOld =NULL;
  p->dx =NULL;
  p->flag2=NULL;
  p->flag3=NULL;
  p->flag4=NULL;
  p->dq =NULL;
  p->lt =NULL;
  p->tq =NULL;
  p->ltq =NULL;
  p->nTime =NULL;
  p->nSD =NULL;
  p->x =NULL;
  p->q =NULL;
  p->diffxq =NULL;
  p->SD =NULL;
  p->ZS =NULL;
  p->nZS =NULL;
  p->nextEventTime = NULL;
  p->zcSign = NULL;
  return (p);
};

void LIQSS2_init (LIQSS2_quantizerState p, const uinteger STATES){
    p->tx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->a = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->mpr = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->u0 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->u1 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lqu = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->qAux = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->oldDx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lquOld = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->dx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->flag2 = (modelica_boolean*)malloc(STATES * sizeof(modelica_boolean));
    p->flag3 = (modelica_boolean*)malloc(STATES * sizeof(modelica_boolean));
    p->flag4 = (modelica_boolean*)malloc(STATES * sizeof(modelica_boolean));
    p->dq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lt = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->tq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->nSD = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
    p->nTime = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->ltq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->nZS = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
    p->nextEventTime = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->zcSign = (modelica_real*)malloc(STATES * sizeof(modelica_real));

    uinteger i=0;
    p->x = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->q = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->diffxq = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->SD = (modelica_integer**) malloc(sizeof(modelica_integer*)*STATES);
    p->ZS = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    for(i=0; i<STATES; i++){
        p->x[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
        p->q[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order));
        p->diffxq[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
        p->SD[i] = (modelica_integer*) malloc(sizeof(modelica_integer)*(STATES));
        p->ZS[i] = (modelica_real*) malloc(sizeof(modelica_real)*(STATES));
    }
}

void LIQSS2_freeQuantizer (LIQSS2_quantizerState p, const uinteger STATES){
	uinteger i=0;
    for(i=0; i<STATES; i++){
        free(p->x[i]);
		free(p->q[i]);
		free(p->diffxq[i]);
		free(p->SD[i]);
        free(p->ZS[i]);
	}
	free(p->a);
	free(p->mpr);
	free(p->u0);
	free(p->u1);
	free(p->lqu);
	free(p->qAux);
	free(p->oldDx);
	free(p->lquOld);
	free(p->dx);
	free(p->flag2);
	free(p->flag3);
	free(p->flag4);
	free(p->dq);
	free(p->lt);
	free(p->tq);
	free(p->nSD);
	free(p->nTime);
	free(p->ltq);
	free(p->tx);
    free(p->nZS);
    free(p->nextEventTime);
}


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
            // printf("disc: %.20f\n", disc);
            if(disc < 0)
				return 1.0 /0.0;
			else{
				modelica_real sd = sqrt(disc);
				modelica_real r1 = (-diffxq[index][1] + sd) / (2 * diffxq[index][2]);
                // printf("r1: %.20f\n", r1);
                if(r1 > 0)
					mpr = r1;
				else
					mpr = 1.0/0.0;

				r1 = (-diffxq[index][1] - sd) / (2 * diffxq[index][2]);
                // printf("r2: %.20f\n", r1);
                if(r1 > 0 && r1 < mpr){
					mpr = r1;
				}
				return mpr;
			}
		}
	}
	else
		exit(1);
}





static modelica_real minRootPosCoeff(modelica_real* diffxq, const uinteger index, const uinteger order){
	modelica_real mpr = 0;
	if(order == 1){
		if(diffxq[1] == 0)
			return 1.0 /0.0;
		else{
			mpr = -diffxq[0] / diffxq[1];
			if(mpr < 0)
				return 1.0 /0.0;
			else
				return mpr;
		}
	}
	else
		exit(1);
}





static modelica_real minRootPosNoIndex(modelica_real* diffxq, const uinteger order){
	modelica_real mpr = 0;
	if(order == 1){
		if(diffxq[1] == 0)
			return 1.0 /0.0;
		else{
			mpr = -diffxq[0] / diffxq[1];
			if(mpr < 0)
				return 1.0 /0.0;
			else
				return mpr;
		}
	}
	else if(order == 2){
		if(diffxq[2] == 0 || (1000*fabs(diffxq[2]) < fabs(diffxq[1]))){
			if(diffxq[1] == 0)
				return 1.0 /0.0;
			else{
				mpr = -diffxq[0] / diffxq[1];
				if(mpr < 0)
					return 1.0 / 0.0;
				else
					return mpr;
			}
		}
		else{
			modelica_real disc = diffxq[1] * diffxq[1] - 4 * diffxq[2] * diffxq[0];
			if(disc < 0)
				return 1.0 /0.0;
			else{
				modelica_real sd = sqrt(disc);
				modelica_real r1 = (-diffxq[1] + sd) / (2 * diffxq[2]);
				if(r1 > 0)
					mpr = r1;
				else
					mpr = 1.0/0.0;
				r1 = (-diffxq[1] - sd) / (2 * diffxq[2]);
				if(r1 > 0 && r1 < mpr){
					mpr = r1;
				}
				return mpr;
			}
		}
	}
	else
		exit(1);
}


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
	return x[index][0] + dt * x[index][1] + dt * dt * x[index][2];
}

static modelica_real calcDerivative(modelica_real** x, const uinteger index, const modelica_real dt){
	return x[index][1] + 2 * dt * x[index][2];
}

static modelica_real calcQState(modelica_real** q, const uinteger index, const modelica_real dt){
	return q[index][0] + dt * q[index][1];
}

static void LIQSS2_recomputeNextTimes(uinteger STATES, LIQSS2_quantizerState p, modelica_real t, modelica_real ft, int dt_index){
    modelica_real diffQ = 0, timeaux = 0;
    uinteger i =0, j =0;
    for(j=0;j<p->nSD[dt_index];j++){
        i= p->SD[dt_index][j];
    // for(j=0;j<STATES;j++){
    //     i=j;
        // if(t>0.00001 && t<0.00001001){
        //     printf("i:  %d\n", i);
        //     printf("t:  %.20f\n", t);
        //     printf("p->nTime[i]:   %.20f \n", p->nTime[i]);
        // }
        if(1){
            // printf("LIQSS2_recomputeNextTimes\n");
            // printf("i:   %d\n", i);
            // printf("t:   %.20f\n", t);
            // printf("p->nTime[i]:   %.20f\n", p->nTime[i]);
            // printf("p->q[i][1]:   %.20f \n", p->q[i][1]);
            // printf("p->x[i][1]:   %.20f \n", p->x[i][1]);
        }
        if(p->ltq[i] == t){
            diffQ = p->q[i][0] - p->qAux[i];
            if(fabs(diffQ) > p->lqu[i] * 1e-6){
                p->a[i] = (p->x[i][1] - p->oldDx[i]) / diffQ;
                if(p->a[i] > 0)
                    p->a[i] = 0;

                // printf("! a[var]: %.20f\n", p->a[i]);
                // printf("! x[cf1]: %.20f\n", p->x[i][1]);
                // printf("! oldDx[var]: %.20f\n", p->oldDx[i]);
                // printf("! diffQ: %.20f\n", diffQ);
            }
        }
        else
            p->flag3[i] = false;

        p->u0[i] = p->x[i][1] - p->q[i][0] * p->a[i];
        p->u1[i] = 2 * p->x[i][2] - p->q[i][1] * p->a[i];
        p->lt[i] = t;
        if(p->flag4[i])
            p->nTime[i] = t;
        else{
            p->diffxq[i][1] = p->q[i][1] - p->x[i][1];
        //     if(1){
            // printf("p->q[i][1] %.20f\n", p->q[i][1]);
            // printf("p->x[i][1] %.20f\n", p->x[i][1]);
        //     printf("p->diffxq[i][1] %.20f\n", p->diffxq[i][1]);
        // }
            p->diffxq[i][2] = -p->x[i][2];
            p->diffxq[i][0] = p->q[i][0] - p->dq[i] + p->lqu[i] - p->x[i][0];
            p->nTime[i] = t + minRootPos(p->diffxq,i,2);
            // printf("nTime[var]: %.20f\n", p->nTime[i]);
            p->diffxq[i][0] = p->q[i][0] - p->dq[i] - p->lqu[i] - p->x[i][0];
            timeaux = t + minRootPos(p->diffxq,i,2);
            // printf("timeaux: %.20f\n", timeaux);
            if(timeaux < p->nTime[i])
                p->nTime[i] = timeaux;
            if(p->a[i] != 0 && fabs(p->x[i][2]) > 1e-10 && !p->flag3[i] && !p->flag2[i]){
                modelica_real coeff[2];
                coeff[0] = p->a[i] * p->a[i] * p->q[i][0] + p->a[i] * p->u0[i] + p->u1[i];
                coeff[1] = p->a[i] * p->a[i] * p->q[i][1] + p->a[i] * p->u1[i];
                // printf("coeff[0] %.20f\n", coeff[0]);
                // printf("coeff[1] %.20f\n", coeff[1]);
                timeaux = t + minRootPosCoeff(&coeff, i, 1);
                // printf("??timeaux: %.20f\n", timeaux);
                if(timeaux < p->nTime[i]){
                    p->flag2[i] = true;
                    p->nTime[i] = timeaux;
                    p->lquOld[i] = p->lqu[i];
                }
            }
            else
                p->flag2[i] = false;

            if(p->nTime[i] > ft)
                p->nTime[i] = ft;

            modelica_real err1 = p->q[i][0] - p->x[i][0] + (p->diffxq[i][1] * (p->nTime[i] - t) / 2) + p->diffxq[i][2] * pow ((p->nTime[i] - t) / 2, 2);
            if (fabs (err1) > 3 * fabs (p->lqu[i])){
                // printf("var: %d\n", i);
                // printf("err1 %.20f\n", err1);

                // printf("t %.20f\n\n", t);
                // printf("lqu[var] %.20f\n", p->lqu[i]);
        	    p->nTime[i] = t + ft * 1.0e-14;
        	}
        }
        // if(t>0.00001 && t<0.00001001){
        //     //printf("Now we are in LIQSS2_recomputeNextTimes\n");
        //     //printf("We are considering the case where i=6\n");
        //     printf("p->nTime[i]:   %.20f \n", p->nTime[i]);
        //     printf("p->diffxq[i][1]:   %.20f \n", p->diffxq[i][1]);
        //     printf("p->diffxq[i][2]:   %.20f \n", p->diffxq[i][2]);
        //     printf("p->diffxq[i][0]:   %.20f \n", p->q[i][0] - p->dq[i] + p->lqu[i] - p->x[i][0]);
        //     printf("p->diffxq[i][0]:   %.20f \n", p->q[i][0] - p->dq[i] - p->lqu[i] - p->x[i][0]);
        // }
        if(1){
            //printf("Now we are in LIQSS2_recomputeNextTimes\n");
            //printf("We are considering the case where i=6\n");
            // printf("p->q[i][0]:   %.20f \n", p->q[i][0]);
            // printf("p->dq[i]:   %.20f \n", p->dq[i]);
            // printf("p->lqu[i]:   %.20f \n", p->lqu[i]);
            // printf("p->x[i][0]:   %.20f \n", p->x[i][0]);
            // //
            // printf("p->nTime[i]:   %.20f \n", p->nTime[i]);
            // printf("p->diffxq[i][1]:   %.20f \n", p->diffxq[i][1]);
            // printf("p->diffxq[i][2]:   %.20f \n", p->diffxq[i][2]);
            // printf("p->diffxq[i][0]:   %.20f \n", p->q[i][0] - p->dq[i] + p->lqu[i] - p->x[i][0]);
            // printf("p->diffxq[i][0]:   %.20f \n\n", p->q[i][0] - p->dq[i] - p->lqu[i] - p->x[i][0]);
        }
    }
}

static void LIQSS2_printModelData(DATA* data){
    uinteger i = 0;

    SIMULATION_DATA *sDataOld = (SIMULATION_DATA*)data->localData[0];
    modelica_real* realVars = sDataOld->realVars;
    printf("\trealVarsData: (%d)\n", data->modelData->nVariablesReal);
    printf("-------------------------------------------\n");
	for(i=0;i< data->modelData->nVariablesReal;i++){
		printf("%s   %f\n",data->modelData->realVarsData[i].info.name, realVars[i]);
	}
    printf("-------------------------------------------\n");

	printf("\tnVariablesInteger: (%d)\n", data->modelData->nVariablesInteger);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nVariablesInteger; i++){
		printf("%s\n", data->modelData->integerVarsData[i].info.name);
	}
    printf("-------------------------------------------\n");

	printf("\tbooleanVarsData: (%d)\n", data->modelData->nVariablesBoolean);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nVariablesBoolean; i++){
		printf("%s\n", data->modelData->booleanVarsData[i].info.name);
	}
    printf("-------------------------------------------\n");

	printf("\tintegerVarsData: (%d)\n", data->modelData->nVariablesInteger);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nVariablesInteger; i++){
		printf("%s\n", data->modelData->integerVarsData[i].info.name);
	}
    printf("-------------------------------------------\n");

	printf("\tstringVarsData: (%d)\n", data->modelData->nVariablesString);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nVariablesString; i++){
		printf("%s\n", data->modelData->stringVarsData[i].info.name);
	}
    printf("-------------------------------------------\n");

	printf("\tnParametersInteger: (%d)\n", data->modelData->nParametersInteger);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nParametersInteger; i++){
		printf("%s\n", data->modelData->integerParameterData[i].info.name);
	}
    printf("-------------------------------------------\n");

    printf("\trealParameterData: (%d)\n", data->modelData->nParametersReal);
    printf("-------------------------------------------\n");
	for(i = 0; i < data->modelData->nParametersReal; i++){
		printf("%s    %.16f\n", data->modelData->realParameterData[i].info.name, data->simulationInfo->realParameter[i]);
	}
    printf("-------------------------------------------\n");

    printf("\tnStateSets: (%d)\n", data->modelData->nStateSets);
}


static void LIQSS_initialiation(LIQSS2_quantizerState p, modelica_real* state, modelica_real* stateDer, const uinteger STATES, modelica_real ft, modelica_real tolerance, modelica_real absTolerance){
    modelica_real mpr1 = 0;
    modelica_real mpr2 = 0;
    uinteger i = 0;
    for(i=0;i<STATES;i++){
        p->tx[i] = 0;
        p->tq[i] = 0;
        p->ltq[i] = 0;
        p->lt[i] = 0;
        p->flag2[i] = false;
        p->flag3[i] = false;
        p->flag4[i] = false;
        p->x[i][0] = state[i];
        p->x[i][1] = stateDer[i];
        p->q[i][0] = state[i];
        p->x[i][2] = 0;
        p->u0[i] = p->x[i][1] - p->q[i][0] * p->a[i];
        p->u1[i] = 2*p->x[i][2] - p->q[i][1] * p->a[i];
        p->diffxq[i][1] = p->q[i][1] - p->x[i][1];
        p->diffxq[i][2] = -p->x[i][2];
        p->lqu[i] = p->x[i][0] * tolerance;
        if(p->lqu[i] < absTolerance)
            p->lqu[i] = absTolerance;
        p->diffxq[i][0] = p->q[i][0] - p->dq[i] + p->lqu[i] - p->x[i][0];
        mpr1 = minRootPos(p->diffxq,i,2);
        p->diffxq[i][0] = p->q[i][0] - p->dq[i] - p->lqu[i] - p->x[i][0];
        mpr2 = minRootPos(p->diffxq,i,2);

        if(mpr1 < mpr2)
            p->mpr[i] = mpr1;
        else
            p->mpr[i] = mpr2;
        if(p->mpr[i] > ft)
            p->mpr[i] = ft;
    }
}


static void LIQSS2_updateQuantizedState(LIQSS2_quantizerState p, modelica_real t, int dt_index){
    modelica_real elapsed = 0;
    p->flag3[dt_index] = false;
    elapsed = t - p->tq[dt_index];
    p->qAux[dt_index] = p->q[dt_index][0] + elapsed * p->q[dt_index][1];
    p->oldDx[dt_index] = p->x[dt_index][1];
    elapsed = t - p->lt[dt_index];
    p->ltq[dt_index] = t;
    p->u0[dt_index] = p->u0[dt_index] + elapsed * p->u1[dt_index];

    if(p->flag2[dt_index]){
        p->lqu[dt_index] = p->lquOld[dt_index];
        p->flag2[dt_index] = false;
        p->q[dt_index][0] = p->qAux[dt_index];
    }
    else
        p->q[dt_index][0] = p->x[dt_index][0];

    if(p->a[dt_index] < -1e-30){
        // printf("Inside the if\n");
        // printf("a[var]:  %.20f\n", p->a[dt_index]);
    	if(p->x[dt_index][2] < 0){
    		p->dx[dt_index] = p->a[dt_index] * p->a[dt_index] * (p->q[dt_index][0] + p->lqu[dt_index]) + p->a[dt_index] * p->u0[dt_index] + p->u1[dt_index];
    		if(p->dx[dt_index] <= 0)
    			p->dq[dt_index] = p->lqu[dt_index];
    		else{
    			p->dq[dt_index] = (-p->u1[dt_index] / p->a[dt_index] / p->a[dt_index]) - (p->u0[dt_index] / p->a[dt_index]) - p->q[dt_index][0];
    			p->flag3[dt_index] = true;
    			if(fabs(p->dq[dt_index]) > p->lqu[dt_index])
    				p->dq[dt_index] = p->lqu[dt_index];
    		}
    	}
    	else{
    		p->dx[dt_index] = p->a[dt_index] * p->a[dt_index] * (p->q[dt_index][0] - p->lqu[dt_index]) + p->a[dt_index] * p->u0[dt_index] + p->u1[dt_index];
    		if(p->dx[dt_index] >= 0)
    			p->dq[dt_index] = -p->lqu[dt_index];
    		else{
    			p->dq[dt_index] = (-p->u1[dt_index] / p->a[dt_index] / p->a[dt_index]) - (p->u0[dt_index] / p->a[dt_index]) - p->q[dt_index][0];
    			p->flag3[dt_index] = true;
    			if(fabs(p->dq[dt_index]) > p->lqu[dt_index]){
    				p->dq[dt_index] = -p->lqu[dt_index];
    			}
    		}
    	}
    	if(p->q[dt_index][1] * p->x[dt_index][1] < 0 && !p->flag2[dt_index] && !p->flag3[dt_index] && !p->flag4[dt_index]){
    		if(p->q[dt_index][1] < 0)
    			p->dq[dt_index] = p->qAux[dt_index] - p->q[dt_index][0] - fabs(p->lquOld[dt_index]) * 0.1;
    		else
    			p->dq[dt_index] = p->qAux[dt_index] - p->q[dt_index][0] + fabs(p->lquOld[dt_index]) * 0.1;
    		p->flag4[dt_index] = true;
    	}
    	else if (p->flag4[dt_index]){
    		p->flag4[dt_index] = false;
    		if(fabs((-p->u1[dt_index] / p->a[dt_index] / p->a[dt_index]) - (p->u0[dt_index] / p->a[dt_index]) - p->q[dt_index][0]) < 3* p->lqu[dt_index]){
    			p->dq[dt_index] = (-p->u1[dt_index] / p->a[dt_index] / p->a[dt_index]) - (p->u0[dt_index] / p->a[dt_index]) - p->q[dt_index][0];
    			p->flag3[dt_index] = true;
    		}
    	}
    }
    else{
        // printf("Inside the else\n");
    	p->flag4[dt_index] = false;
    	if(p->x[dt_index][2] < 0){
            // printf("Inside the 1st if\n");
    		p->dq[dt_index] = -p->lqu[dt_index];
        }
    	else{
            // printf("Inside the 1st else\n");
    		p->dq[dt_index] = p->lqu[dt_index];
        }
    }

    if(fabs(p->dq[dt_index]) > 2* p->lqu[dt_index]){
        // printf("Inside the second if\n");
        if(p->dq[dt_index] > 0){
            // printf("Inside the 2nd if\n");
            p->dq[dt_index] = p->lqu[dt_index];
        }
        else{
            // printf("Inside the 2nd else\n");
            p->dq[dt_index] = -p->lqu[dt_index];
        }
    }

    p->q[dt_index][0] = p->q[dt_index][0] + p->dq[dt_index];
    if(p->flag3[dt_index]){
        // printf("Inside the third if\n");
        p->q[dt_index][1] = p->a[dt_index] * p->q[dt_index][0] + p->u0[dt_index];
    }
    else{
        // printf("Inside the third else\n");
        p->q[dt_index][1] = p->x[dt_index][1];
    }

    // printf("Inside update\n");
    // printf("dt_index: %d\n", dt_index);
    // printf("p->dq[dt_index]:   %.20f\n", p->dq[dt_index]);
    // printf("p->lqu[dt_index]:   %.20f\n\n", p->lqu[dt_index]);
}


static int LIQSS_write(LIQSS2_quantizerState quantizer, DATA* data, int type, modelica_real t, const uinteger STATES, int dt_index){
    switch(type){
        case 1:
        {
            uinteger i = 0;
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
            return 0;
        }
        break;
        case 2:
        {
            uinteger i = 0;
            char buf[100];
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
            return 0;
        }
        break;
        case 3:
        {
            char buf[100];
            sprintf(buf, "%s%s.csv", "/home/anton/Desktop/Results/", data->modelData->realVarsData[dt_index].info.name);
			FILE *fs = fopen(buf, "a");
			if(fs == NULL){
				printf("Couldn't open file\n");
				return 1;
			}
			fprintf(fs,"%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",t,quantizer->q[dt_index][0],quantizer->q[dt_index][1],quantizer->q[dt_index][2],quantizer->x[dt_index][0],quantizer->x[dt_index][1],quantizer->x[dt_index][2]);
			fclose(fs);
            return 0;
        }
        break;
        default:
            printf("Invalid type.\n");
            return 1;
    }
}


static void ZeroCrossing(LIQSS2_quantizerState quantizer, modelica_real *coeff, int ZSIndex){
    coeff[0] = quantizer->q[ZSIndex][0] - 10.0;
    coeff[1] = quantizer->q[ZSIndex][1];
}

static int signQSS(modelica_real f){
    if(f >= 0)
        return 1;
    else
        return (-1);
}


static void FRW_nextEventTime(LIQSS2_quantizerState p, modelica_real t, int dt_index){
    int i, j, s;
    modelica_real e = 0;
    int order = p->order - 1;
    modelica_real *coeff = NULL;
    coeff = (modelica_real*)calloc(2, sizeof(modelica_real));


    for (i = 0; i < p->nZS[dt_index]; i++){
        j = p->ZS[dt_index][i];
        e = t - p->tq[dt_index];
        if(e>0){
            // quantizer->q[ZSIndex][0] = calcState(quantizer->q, ZSIndex, e);
    		// quantizer->q[ZSIndex][1] = calcDerivative(quantizer->q, ZSIndex, e);
            printf("Hierdie moenie eintlik nou gebeur nie\n");
        }
    }
    ZeroCrossing(p, coeff, dt_index);
    s = signQSS(coeff[0]);

    // Vet hack hier wat volg:
    p->zcSign[dt_index] = s;

    modelica_real zcHyst = 1e-10;
    coeff[0] += p->zcSign[dt_index] + zcHyst;
    p->nextEventTime[dt_index] = t + minRootPosNoIndex(coeff, 1);
}
