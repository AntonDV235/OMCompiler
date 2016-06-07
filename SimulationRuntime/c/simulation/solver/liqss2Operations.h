void * checkedMalloc (unsigned long long len);

LIQSS2_quantizerState LIQSS2_QuantizerState ();

void LIQSS2_init (LIQSS2_quantizerState p, const uinteger STATES);

void LIQSS2_freeQuantizer (LIQSS2_quantizerState p, const uinteger STATES);

static modelica_real minRootPos(modelica_real** diffxq, const uinteger index, const uinteger order);

static modelica_real calcState(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcDerivative(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcQState(modelica_real** q, const uinteger index, const modelica_real dt);
// static modelica_boolean calcAllDerivatives(modelica_real** q, modelica_real** x, DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo);
static modelica_boolean LIQSS2_calculateState(DATA* data, threadData_t *threadData);

static void LIQSS2_recomputeNextTimes(LIQSS2_quantizerState p, modelica_real t, modelica_real ft, int dt_index);

static void LIQSS2_updateQuantizedState(LIQSS2_quantizerState p, modelica_real t, int dt_index);
