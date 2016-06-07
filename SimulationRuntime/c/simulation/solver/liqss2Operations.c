#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solver_main.h"
#include "simulation/simulation_runtime.h"
#include "simulation/results/simulation_result.h"
#include "openmodelica_func.h"
#include "util/omc_error.h"
#include "simulation/options.h"
// #include "simulation/solver/liqss2Operations.h"
// #include "simulation/solver/liqss2.h"

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
  printf("p->order: %d\n", p->order);
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
  return (p);
};

void LIQSS2_init (LIQSS2_quantizerState p, const uinteger STATES){
    p->tx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    //fail = (tx == NULL) ? 1 : ( 0 | fail);
    p->a = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->mpr = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->u0 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->u1 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lqu = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->qAux = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->oldDx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lquOld = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->dx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->flag2 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
    p->flag3 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
    p->flag4 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
    p->dq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->lt = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->tq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->nSD = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
    p->nTime = (modelica_real*)malloc(STATES * sizeof(modelica_real));
    p->ltq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
	printf("p->order: %d\n", p->order);
    p->order=5;
	printf("p->order: %d\n", p->order);

    uinteger i=0;
    p->x = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->q = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->diffxq = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
    p->SD = (modelica_integer**) malloc(sizeof(modelica_integer*)*STATES);
    for(i=0; i<STATES; i++){
        p->x[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
        p->q[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order));
        p->diffxq[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
        p->SD[i] = (modelica_integer*) malloc(sizeof(modelica_integer)*(STATES));
    }
}

void LIQSS2_freeQuantizer (LIQSS2_quantizerState p, const uinteger STATES){
	uinteger i=0;
    for(i=0; i<STATES; i++){
        free(p->x[i]);
		free(p->q[i]);
		free(p->diffxq[i]);
		free(p->SD[i]);
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
}



//
// LIQSS2_quantizerState LIQSS2_QuantizerState (){
//   LIQSS2_quantizerState p = checkedMalloc (sizeof(*p));
//   p->order = 2;
//   printf("p->order: %d\n", p->order);
//   p->x = NULL;
//   p->diffxq = NULL;
//   p->q = NULL;
//   p->SD =NULL;
//   p->tx = NULL;
//   p->a = NULL;
//   p->mpr =NULL;
//   p->u0 =NULL;
//   p->u1 =NULL;
//   p->lqu =NULL;
//   p->qAux =NULL;
//   p->oldDx =NULL;
//   p->lquOld =NULL;
//   p->dx =NULL;
//   p->flag2=NULL;
//   p->flag3=NULL;
//   p->flag4=NULL;
//   p->dq =NULL;
//   p->lt =NULL;
//   p->tq =NULL;
//   p->ltq =NULL;
//   p->nTime =NULL;
//   p->nSD =NULL;
//   p->x =NULL;
//   p->q =NULL;
//   p->diffxq =NULL;
//   p->SD =NULL;
//   return (p);
// };
//
// void LIQSS2_init (LIQSS2_quantizerState p, const uinteger STATES){
//     p->tx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     //fail = (tx == NULL) ? 1 : ( 0 | fail);
//     p->a = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->mpr = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->u0 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->u1 = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->lqu = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->qAux = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->oldDx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->lquOld = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->dx = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->flag2 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
//     p->flag3 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
//     p->flag4 = (modelica_real*)malloc(STATES * sizeof(modelica_boolean));
//     p->dq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->lt = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->tq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->nSD = (modelica_integer*)malloc(STATES * sizeof(modelica_integer));
//     p->nTime = (modelica_real*)malloc(STATES * sizeof(modelica_real));
//     p->ltq = (modelica_real*)malloc(STATES * sizeof(modelica_real));
// 	printf("p->order: %d\n", p->order);
//     p->order=5;
// 	printf("p->order: %d\n", p->order);
//
//     uinteger i=0;
//     p->x = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
//     p->q = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
//     p->diffxq = (modelica_real**) malloc(sizeof(modelica_real*)*STATES);
//     p->SD = (modelica_integer**) malloc(sizeof(modelica_integer*)*STATES);
//     for(i=0; i<STATES; i++){
//         p->x[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
//         p->q[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order));
//         p->diffxq[i] = (modelica_real*) malloc(sizeof(modelica_real)*(p->order+1));
//         p->SD[i] = (modelica_integer*) malloc(sizeof(modelica_integer)*(STATES));
//     }
// }
//
// void LIQSS2_freeQuantizer (LIQSS2_quantizerState p, const uinteger STATES){
// 	uinteger i=0;
//     for(i=0; i<STATES; i++){
//         free(p->x[i]);
// 		free(p->q[i]);
// 		free(p->diffxq[i]);
// 		free(p->SD[i]);
// 	}
// 	free(p->a);
// 	free(p->mpr);
// 	free(p->u0);
// 	free(p->u1);
// 	free(p->lqu);
// 	free(p->qAux);
// 	free(p->oldDx);
// 	free(p->lquOld);
// 	free(p->dx);
// 	free(p->flag2);
// 	free(p->flag3);
// 	free(p->flag4);
// 	free(p->dq);
// 	free(p->lt);
// 	free(p->tq);
// 	free(p->nSD);
// 	free(p->nTime);
// 	free(p->ltq);
// 	free(p->tx);
// }
