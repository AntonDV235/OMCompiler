// static modelica_real ddx(modelica_real** x, const uinteger index);
static modelica_real calcState(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcDerivative(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcQState(modelica_real** q, const uinteger index, const modelica_real dt);
static modelica_boolean calcAllDerivatives(modelica_real** q, modelica_real** x, DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo);
static modelica_real minRootPos(modelica_real** diffxq, const uinteger index, const uinteger order);
static modelica_boolean LIQSS2_calculateState(DATA* data, threadData_t *threadData);
static int prefixedName_LIQSS2Simulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo);


typedef struct LIQSS2_quantizerState_ *LIQSS2_quantizerState;

struct LIQSS2_quantizerState_
{
  int order;
  modelica_real** x; //!<
  modelica_real** diffxq; //!<
  modelica_real** q; //!<
  modelica_integer** SD;
  modelica_real* tx;
  modelica_real* a;
  modelica_real* mpr;
  modelica_real* u0;
  modelica_real* u1;
  modelica_real* lqu;
  modelica_real* qAux;
  modelica_real* oldDx;
  modelica_real* lquOld;
  modelica_real* dx;
  modelica_boolean* flag2;
  modelica_boolean* flag3;
  modelica_boolean* flag4;
  modelica_real* dq;
  modelica_real* lt;
  modelica_real* tq;
  modelica_real* ltq;
  modelica_real* nTime;
  modelica_integer* nSD;
};
