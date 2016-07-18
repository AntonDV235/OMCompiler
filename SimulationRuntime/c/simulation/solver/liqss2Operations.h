void * checkedMalloc (unsigned long long len);

LIQSS2_quantizerState LIQSS2_QuantizerState ();

void LIQSS2_init (LIQSS2_quantizerState p, const uinteger STATES);

void LIQSS2_freeQuantizer (LIQSS2_quantizerState p, const uinteger STATES);

static modelica_real minRootPos(modelica_real** diffxq, const uinteger index, const uinteger order);

// static modelica_real minRootPosCoeff(modelica_real* diffxq, const uinteger index, const uinteger order);

// static modelica_real minRootPosNoIndex(modelica_real* diffxq, const uinteger order);

static modelica_real calcState(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcDerivative(modelica_real** x, const uinteger index, const modelica_real dt);
static modelica_real calcQState(modelica_real** q, const uinteger index, const modelica_real dt);
// static modelica_boolean calcAllDerivatives(modelica_real** q, modelica_real** x, DATA* data, threadData_t *threadData, SOLVER_INFO* solverInfo);
static modelica_boolean LIQSS2_calculateState(DATA* data, threadData_t *threadData);

static void LIQSS2_recomputeNextTimes(uinteger STATES, LIQSS2_quantizerState p, modelica_real t, modelica_real ft, int dt_index);

static void LIQSS2_updateQuantizedState(LIQSS2_quantizerState p, modelica_real t, int dt_index);

static void LIQSS_initialiation(LIQSS2_quantizerState p, modelica_real* state, modelica_real* stateDer, const uinteger STATES, modelica_real ft, modelica_real tolerance, modelica_real absTolerance);

static void LIQSS2_printModelData(DATA* data);

// Type 1 -- Creating StateVars.txt
// Type 2 -- For writing column headers
// Type 3 -- For writing simulation results
static int LIQSS_write(LIQSS2_quantizerState quantizer, DATA* data, int type, modelica_real t, const uinteger STATES, int dt_index);

static void FRW_nextEventTime(LIQSS2_quantizerState quantizer, modelica_real t, int dt_index);

static void ZeroCrossing(LIQSS2_quantizerState quantizer, modelica_real *coeff, int ZSIndex);

static int signQSS(modelica_real f);
