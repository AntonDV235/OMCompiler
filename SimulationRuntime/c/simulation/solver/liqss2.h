static modelica_real ddx(DATA* data,modelica_real t, modelica_real** x, const uinteger index, modelica_real* state);

static modelica_real dx(DATA* data, modelica_real t, modelica_real** x, const uinteger index, modelica_real* state);

static int prefixedName_LIQSS2Simulation(DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo);

typedef struct LIQSS2_quantizerState_ *LIQSS2_quantizerState;

struct LIQSS2_quantizerState_
{
  int order;
  int events; //!<
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
  modelica_integer* nZS;
  modelica_real** ZS;
  modelica_real* nextEventTime;
  modelica_real* zcSign;
};
