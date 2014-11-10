#include <Core/Modelica.h>
#include <Solver/Peer/Peer.h>
#include <Core/Math/Functions.h>

Peer::Peer(IMixedSystem* system, ISolverSettings* settings)
    : SolverDefaultImplementation(system, settings),
      _peersettings(dynamic_cast<ISolverSettings*>(_settings))
/*      _cvodeMem(NULL),
      _z(NULL),
      _zInit(NULL),
      _zWrite(NULL),
      _dimSys(0),
      _cv_rt(0),
      _outStps(0),
      _locStps(0),
      _idid(0),
      _hOut(0.0),
      _tOut(0.0),
      _tZero(0.0),
      _zeroSign(NULL),
      _absTol(NULL),
      _cvode_initialized(false),
      _tLastEvent(0.0),
      _event_n(0),
      _properties(NULL),
      _continuous_system(NULL),
      _event_system(NULL),
      _mixed_system(NULL),
      _time_system(NULL),
    _delta(NULL),
    _ysave(NULL) */
{
  //_data = ((void*) this);
}

Peer::~Peer()
{
  /*
  if (_z)
    delete[] _z;
  if (_zInit)
    delete[] _zInit;
  if (_zeroSign)
    delete[] _zeroSign;
  if (_absTol)
    delete[] _absTol;
  if (_zWrite)
      delete[] _zWrite;
  if (_cvode_initialized)
  {
    N_VDestroy_Serial(_CV_y0);
    N_VDestroy_Serial(_CV_y);
    N_VDestroy_Serial(_CV_yWrite);
    N_VDestroy_Serial(_CV_absTol);
    CVodeFree(&_cvodeMem);
  }


  if (_sparsePattern_index)
    delete [] _sparsePattern_leadindex;
  if (_sparsePattern_colorCols)
    delete [] _sparsePattern_colorCols;
  if(_sparsePattern_index)
    delete [] _sparsePattern_index;
  if(_delta)
    delete [] _delta;
  if(_ysave)
    delete [] _ysave; */

}

void Peer::initialize()
{
  /*
  _properties = dynamic_cast<ISystemProperties*>(_system);
  _continuous_system = dynamic_cast<IContinuous*>(_system);
  _event_system = dynamic_cast<IEvent*>(_system);
  _mixed_system = dynamic_cast<IMixedSystem*>(_system);
  _time_system = dynamic_cast<ITime*>(_system);
  IGlobalSettings* global_settings = dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings();
  // Kennzeichnung, dass initialize()() (vor der Integration) aufgerufen wurde
  _idid = 5000;
  _tLastEvent = 0.0;
  _event_n = 0;
  SolverDefaultImplementation::initialize();
  _dimSys = _continuous_system->getDimContinuousStates();
  _dimZeroFunc = _event_system->getDimZeroFunc();

  if (_dimSys <= 0)
  {
    _idid = -1;
    throw std::invalid_argument("Peer::initialize()");
  }
  else
  {
    // Allocate state vectors, stages and temporary arrays
    if (_z)
      delete[] _z;
    if (_zInit)
      delete[] _zInit;
    if (_zWrite)
      delete[] _zWrite;
    if (_zeroSign)
      delete[] _zeroSign;
    if (_absTol)
      delete[] _absTol;
  if(_delta)
    delete [] _delta;
    if(_ysave)
    delete [] _ysave;

    _z = new double[_dimSys];
    _zInit = new double[_dimSys];
    _zWrite = new double[_dimSys];
    _zeroSign = new int[_dimZeroFunc];
    _absTol = new double[_dimSys];
  _delta =new double[_dimSys];
  _ysave =new double[_dimSys];

    memset(_z, 0, _dimSys * sizeof(double));
    memset(_zInit, 0, _dimSys * sizeof(double));
  memset(_ysave, 0, _dimSys * sizeof(double));

    // Counter initialisieren
    _outStps = 0;

    if (_cvodesettings->getDenseOutput())
    {
      // Ausgabeschrittweite
      _hOut = global_settings->gethOutput();

    }

    // Allocate memory for the solver
    _cvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void*) _cvodeMem, "CVodeCreate", 0))
    {
      _idid = -5;
      throw std::invalid_argument("Peer::initialize()");
    }

    //
    // Make Peer ready for integration
    //

    // Set initial values for CVODE
    _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
    _continuous_system->getContinuousStates(_zInit);
    memcpy(_z, _zInit, _dimSys * sizeof(double));

    // Get nominal values
    _continuous_system->getNominalStates(_absTol);
    for (int i = 0; i < _dimSys; i++)
      _absTol[i] *= dynamic_cast<ISolverSettings*>(_cvodesettings)->getATol();

    _CV_y0 = N_VMake_Serial(_dimSys, _zInit);
    _CV_y = N_VMake_Serial(_dimSys, _z);
    _CV_yWrite = N_VMake_Serial(_dimSys, _zWrite);
    _CV_absTol = N_VMake_Serial(_dimSys, _absTol);

    if (check_flag((void*) _CV_y0, "N_VMake_Serial", 0))
    {
      _idid = -5;
      throw std::invalid_argument("Peer::initialize()");
    }

    // Initialize Peer (Initial values are required)
    _idid = CVodeInit(_cvodeMem, CV_fCallback, _tCurrent, _CV_y0);
    if (_idid < 0)
    {
      _idid = -5;
      throw std::invalid_argument("Peer::initialize()");
    }

    // Set Tolerances
    _idid = CVodeSVtolerances(_cvodeMem, dynamic_cast<ISolverSettings*>(_cvodesettings)->getRTol(), _CV_absTol);    // RTOL and ATOL
    if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");

    // Set the pointer to user-defined data
    _idid = CVodeSetUserData(_cvodeMem, _data);
    if (_idid < 0)
      throw std::invalid_argument("Peer::initialize()");

    _idid = CVodeSetInitStep(_cvodeMem, 1e-6);    // INITIAL STEPSIZE
    if (_idid < 0)
      throw std::invalid_argument("Peer::initialize()");

    _idid = CVodeSetMaxOrd(_cvodeMem, 5);       // Max Order
    if (_idid < 0)
      throw std::invalid_argument("CVoder::initialize()");

    _idid = CVodeSetMaxConvFails(_cvodeMem, 100);       // Maximale Fehler im Konvergenztest
    if (_idid < 0)
      throw std::invalid_argument("CVoder::initialize()");

    _idid = CVodeSetStabLimDet(_cvodeMem, TRUE);       // Stability Detection
    if (_idid < 0)
      throw std::invalid_argument("CVoder::initialize()");

    _idid = CVodeSetMinStep(_cvodeMem, dynamic_cast<ISolverSettings*>(_cvodesettings)->getLowerLimit());       // MINIMUM STEPSIZE
    if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");

    _idid = CVodeSetMaxStep(_cvodeMem, global_settings->getEndTime() / 10.0);       // MAXIMUM STEPSIZE
    if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");

    _idid = CVodeSetMaxNonlinIters(_cvodeMem, 5);      // Max number of iterations
    if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");
    _idid = CVodeSetMaxErrTestFails(_cvodeMem, 100);
    if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");

    _idid = CVodeSetMaxNumSteps(_cvodeMem, 1e3);            // Max Number of steps
    if (_idid < 0)
      throw std::invalid_argument("Peer::initialize()");

    // Initialize linear solver
  _idid = CVDense(_cvodeMem, _dimSys);
    if (_idid < 0)
      throw std::invalid_argument("Peer::initialize()");

  // Use own jacobian matrix
  //_idid = CVDlsSetDenseJacFn(_cvodeMem, &CV_JCallback);
  if (_idid < 0)
      throw std::invalid_argument("CVode::initialize()");

    if (_dimZeroFunc)
    {
      _idid = CVodeRootInit(_cvodeMem, _dimZeroFunc, &CV_ZerofCallback);

      memset(_zeroSign, 0, _dimZeroFunc * sizeof(int));
      _idid = CVodeSetRootDirection(_cvodeMem, _zeroSign);
      if (_idid < 0)
        throw std::invalid_argument("CVode::initialize()");
      memset(_zeroSign, -1, _dimZeroFunc * sizeof(int));
      memset(_zeroVal, -1, _dimZeroFunc * sizeof(int));

    }

    initializeColoredJac();
    _cvode_initialized = true;

    //
    // CVODE is ready for integration
    //
    // BOOST_LOG_SEV(cvode_lg::get(), cvode_info) << "CVode initialized";
  }
  */
}

void Peer::solve(const SOLVERCALL action)
{
  std::cerr << "using Peer" << std::endl;
  /*
  bool writeEventOutput = (_settings->getGlobalSettings()->getOutputPointType() == ALL);
  bool writeOutput = !(_settings->getGlobalSettings()->getOutputPointType() == EMPTY2);

  if (_cvodesettings && _system)
  {
    // Solver und System fÃ¼r Integration vorbereiten
    if ((action & RECORDCALL) && (action & FIRST_CALL))
    {
      initialize();
      if (writeOutput)
        writeToFile(0, _tCurrent, _h);
      _tLastWrite = 0;

      return;
    }

    if ((action & RECORDCALL) && !(action & FIRST_CALL))
    {
      writeToFile(_accStps, _tCurrent, _h);
      return;
    }

    // Nach einem TimeEvent wird der neue Zustand recorded
    if (action & RECALL)
    {
      _firstStep = true;
      if (writeEventOutput)
        writeToFile(0, _tCurrent, _h);
      if (writeOutput)
        writeCVodeOutput(_tCurrent, _h, _locStps);
    }

    // Solver soll fortfahren
    _solverStatus = ISolver::CONTINUE;

    while (_solverStatus & ISolver::CONTINUE)
    {
      // Zuvor wurde initialize aufgerufen und hat funktioniert => RESET IDID
      if (_idid == 5000)
        _idid = 0;

      // Solveraufruf
      if (_idid == 0)
      {
        // ZÃ¤hler zurÃ¼cksetzen
        _accStps = 0;
        _locStps = 0;

        // Solverstart
        CVodeCore();

      }

      // Integration war nicht erfolgreich und wurde auch nicht vom User unterbrochen
      if (_idid != 0 && _idid != 1)
      {
        _solverStatus = ISolver::SOLVERERROR;
        //throw std::invalid_argument(_idid,_tCurrent,"CVode::solve()");
        throw std::invalid_argument("CVode::solve()");
      }

      // Abbruchkriterium (erreichen der Endzeit)
      else if ((_tEnd - _tCurrent) <= dynamic_cast<ISolverSettings*>(_cvodesettings)->getEndTimeTol())
        _solverStatus = DONE;
    }

    _firstCall = false;

  }
  else
  {

    throw std::invalid_argument("CVode::solve()");
  }
  */
}

void Peer::PeerCore()
{
  /*
  _idid = CVodeReInit(_cvodeMem, _tCurrent, _CV_y);
  _idid = CVodeSetStopTime(_cvodeMem, _tEnd);
  _idid = CVodeSetInitStep(_cvodeMem, 1e-12);
  if (_idid < 0)
    throw std::runtime_error("CVode::ReInit");

  bool writeEventOutput = (_settings->getGlobalSettings()->getOutputPointType() == ALL);
  bool writeOutput = !(_settings->getGlobalSettings()->getOutputPointType() == EMPTY2);

  while (_solverStatus & ISolver::CONTINUE)
  {
    _cv_rt = CVode(_cvodeMem, _tEnd, _CV_y, &_tCurrent, CV_ONE_STEP);

    _idid = CVodeGetNumSteps(_cvodeMem, &_locStps);
    if (_idid != CV_SUCCESS)
      throw std::runtime_error("CVodeGetNumSteps failed. The cvode mem pointer is NULL");

    _idid = CVodeGetLastStep(_cvodeMem, &_h);
    if (_idid != CV_SUCCESS)
      throw std::runtime_error("CVodeGetLastStep failed. The cvode mem pointer is NULL");

    //Check if there was at least one output-point within the last solver interval
    //  -> Write output if true
    if (writeOutput)
      writeCVodeOutput(_tCurrent, _h, _locStps);

    //set completed step to system and check if terminate was called
    if(_continuous_system->stepCompleted(_tCurrent))
        _solverStatus = DONE;

    // Perform state selection
    bool state_selection = stateSelection();
    if (state_selection)
      _continuous_system->getContinuousStates(_z);

    _zeroFound = false;

    // Check if step was successful
    if (check_flag(&_cv_rt, "CVode", 1))
    {
      _solverStatus = ISolver::SOLVERERROR;
      break;
    }

    // A root was found
    if ((_cv_rt == CV_ROOT_RETURN))
    {
      // CVode is setting _tCurrent to the time where the first event occurred
      double _abs = fabs(_tLastEvent - _tCurrent);
      _zeroFound = true;

      if ((_abs < 1e-3) && _event_n == 0)
      {
        _tLastEvent = _tCurrent;
        _event_n++;
      }
      else if ((_abs < 1e-3) && (_event_n >= 1 && _event_n < 500))
      {
        _event_n++;
      }
      else if ((_abs >= 1e-3))
      {
        //restart event counter
        _tLastEvent = _tCurrent;
        _event_n = 0;
      }
      else
        throw std::runtime_error("Number of events exceeded  in time interval " + boost::lexical_cast<string>(_abs) + " at time " + boost::lexical_cast<string>(_tCurrent));

      // CVode has interpolated the states at time 'tCurrent'
      _time_system->setTime(_tCurrent);

      // To get steep steps in the result file, two value points (P1 and P2) must be added
      //
      // Y |   (P2) X...........
      //   |        :
      //   |        :
      //   |........X (P1)
      //   |---------------------------------->
      //   |        ^                         t
      //        _tCurrent

      // Write the values of (P1)
      if (writeEventOutput)
      {
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        writeToFile(0, _tCurrent, _h);
      }

      _idid = CVodeGetRootInfo(_cvodeMem, _zeroSign);

      for (int i = 0; i < _dimZeroFunc; i++)
        _events[i] = bool(_zeroSign[i]);

      if (_mixed_system->handleSystemEvents(_events))
      {
        // State variables were reinitialized, thus we have to give these values to the cvode-solver
        // Take care about the memory regions, _z is the same like _CV_y
        _continuous_system->getContinuousStates(_z);
      }
    }

    if (_zeroFound || state_selection)
    {
      // Write the values of (P2)
      if (writeEventOutput)
      {
        // If we want to write the event-results, we should evaluate the whole system again
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        writeToFile(0, _tCurrent, _h);
      }

      _idid = CVodeReInit(_cvodeMem, _tCurrent, _CV_y);
      if (_idid < 0)
        throw std::runtime_error("CVode::ReInit()");

      // Der Eventzeitpunkt kann auf der Endzeit liegen (Time-Events). In diesem Fall wird der Solver beendet, da CVode sonst eine interne Warnung schmeißt
      if (_tCurrent == _tEnd)
        _cv_rt = CV_TSTOP_RETURN;
    }

    // ZÃ¤hler fÃ¼r die Anzahl der ausgegebenen Schritte erhÃ¶hen
    ++_outStps;
    _tLastSuccess = _tCurrent;

    if (_cv_rt == CV_TSTOP_RETURN)
    {
      _time_system->setTime(_tEnd);
//      _continuous_system->setContinuousStates(NV_DATA_S(_CV_y));
//      _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
      if(writeOutput)
        writeToFile(0, _tEnd, _h);
      _accStps += _locStps;
      _solverStatus = DONE;
    }
  }
  */
}

void Peer::writePeerOutput(const double &time, const double &h, const int &stp)
{
  /*
  if (stp > 0)
  {
    if (_cvodesettings->getDenseOutput())
    {
      _bWritten = false;
      double *oldValues = NULL;

      //We have to find all output-points within the last solver step
      while (_tLastWrite + dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings()->gethOutput() <= time)
      {
        if (!_bWritten)
        {
          //Rescue the calculated derivatives
          oldValues = new double[_continuous_system->getDimRHS()];
          _continuous_system->getRHS(oldValues);
        }
        _bWritten = true;
        _tLastWrite = _tLastWrite + dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings()->gethOutput();
        //Get the state vars at the output-point (interpolated)
        _idid = CVodeGetDky(_cvodeMem, _tLastWrite, 0, _CV_yWrite);
        _time_system->setTime(_tLastWrite);
        _continuous_system->setContinuousStates(NV_DATA_S(_CV_yWrite));
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        SolverDefaultImplementation::writeToFile(stp, _tLastWrite, h);
      }      //end if time -_tLastWritten
      if (_bWritten)
      {
        _time_system->setTime(time);
        _continuous_system->setContinuousStates(_z);
        _continuous_system->setRHS(oldValues);
        delete[] oldValues;
        //_continuous_system->evaluateAll(IContinuous::CONTINUOUS);
      }
      else if (time == _tEnd && _tLastWrite != time)
      {
        _idid = CVodeGetDky(_cvodeMem, time, 0, _CV_y);
        _time_system->setTime(time);
        _continuous_system->setContinuousStates(NV_DATA_S(_CV_y));
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        SolverDefaultImplementation::writeToFile(stp, _tEnd, h);
      }
    }
    else
      SolverDefaultImplementation::writeToFile(stp, time, h);
  }
  */
}

bool Peer::stateSelection()
{
  return SolverDefaultImplementation::stateSelection();
}
int Peer::calcFunction(const double& time, const double* y, double* f)
{
  /*
  MEASURETIME_REGION_DEFINE(cvodeCalcFunctionHandler, "CVodeCalcFunction");

  if(MeasureTime::getInstance() != NULL)
  {
      MEASURETIME_START(measuredFunctionStartValues, cvodeCalcFunctionHandler, "CVodeCalcFunction");
  }
  int returnValue = 0;
  try
  {
    _time_system->setTime(time);
    _continuous_system->setContinuousStates(y);
    _continuous_system->evaluateODE(IContinuous::CONTINUOUS);
    _continuous_system->getRHS(f);
  }      //workaround until exception can be catch from c- libraries
  catch (std::exception& ex)
  {
    std::string error = ex.what();
    cerr << "CVode integration error: " << error;
    returnValue = 1;
  }

  if(MeasureTime::getInstance() != NULL)
  {
      MEASURETIME_END(measuredFunctionStartValues, measuredFunctionEndValues, measureTimeFunctionsArray[0], cvodeCalcFunctionHandler);
  }
  return returnValue;
  */
  return 0;
}

/*int Peer::calcJacobian(double t, long int N, N_Vector fHelp, N_Vector errorWeight, N_Vector jthCol, double* y, N_Vector fy, DlsMat Jac)
{
  try
  {
    int j,k,l;
  double fnorm, minInc, *f_data, *fHelp_data, *errorWeight_data, h, srur, delta_inv;

  f_data = NV_DATA_S(fy);
  errorWeight_data = NV_DATA_S(errorWeight);
  fHelp_data = NV_DATA_S(fHelp);


  //Get relevant info
  _idid = CVodeGetErrWeights(_cvodeMem, errorWeight);
  if (_idid < 0)
    {
      _idid = -5;
      throw std::invalid_argument("Peer::calcJacobian()");
  }
  _idid = CVodeGetCurrentStep(_cvodeMem, &h);
  if (_idid < 0)
    {
      _idid = -5;
      throw std::invalid_argument("Peer::calcJacobian()");
  }

  srur = sqrt(UROUND);

  fnorm = N_VWrmsNorm(fy, errorWeight);
  minInc = (fnorm != 0.0) ?
           (1000.0 * abs(h) * UROUND * N * fnorm) : 1.0;

  for(j=0;j<N;j++)
  {
    _delta[j] = max(srur*abs(y[j]), minInc/errorWeight_data[j]);
  }

  // Calculation of the jacobian
  for(int i=0; i < _sparsePattern_maxColors; i++)
  {
    for(int ii=0; ii < _dimSys; ii++)
    {
      if((_sparsePattern_colorCols[ii] - 1) == i)
      {
        _ysave[ii] = y[ii];
        y[ii]+= _delta[ii];
      }

    }

    calcFunction(t, y, fHelp_data);

    for(int ii = 0; ii < _dimSys; ii++)
    {
      if((_sparsePattern_colorCols[ii] - 1) == i)
      {

        y[ii] = _ysave[ii];

        if(ii==0)
        {
          j = 0;
        }
        else
        {
          j = _sparsePattern_leadindex[ii-1];

        }
        while(j <_sparsePattern_leadindex[ii])
        {
          l = _sparsePattern_index[j];
          k = l + ii * _dimSys;
          //Jac->data[k] = (fHelp_data[l] - f_data[l])/_delta[l];
          delta_inv = 1.0/_delta[ii];
          Jac->data[k] = (fHelp_data[l] - f_data[l])*delta_inv;
          j++;
        }
      }
    }
  }
  */
  /*
  //Calculation of J without colouring
   for (j = 0; j < N; j++)
   {


    //N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

     _ysave[j] = y[j];

    y[j] += _delta[j];

    calcFunction(t, y, fHelp_data);

    y[j] = _ysave[j];

    delta_inv = 1.0/_delta[j];
    N_VLinearSum(delta_inv, fHelp, -delta_inv, fy, jthCol);

    for(int i=0; i<_dimSys; ++i)
        {
            Jac->data[i+j*_dimSys] = NV_Ith_S(jthCol,i);
        }

    //DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }
  */
/*
 }      //workaround until exception can be catch from c- libraries
  catch (std::exception& ex)
  {
    std::string error = ex.what();
    cerr << "CVode integration error: " << error;
    return 1;
  }


  return 0;
}*/

void Peer::initializeColoredJac()
{
  /*
  _sizeof_sparsePattern_colorCols = _system->getA_sizeof_sparsePattern_colorCols();
  _sparsePattern_colorCols = new int[_sizeof_sparsePattern_colorCols];
  _system->getA_sparsePattern_colorCols( _sparsePattern_colorCols, _sizeof_sparsePattern_colorCols);

  _sizeof_sparsePattern_leadindex = _system->getA_sizeof_sparsePattern_leadindex();
  _sparsePattern_leadindex = new int[_sizeof_sparsePattern_leadindex];
  _system->getA_sparsePattern_leadindex( _sparsePattern_leadindex, _sizeof_sparsePattern_leadindex);


  _sizeof_sparsePattern_index = _system->getA_sizeof_sparsePattern_index();
  _sparsePattern_index = new int[_sizeof_sparsePattern_index];
  _system->getA_sparsePattern_index( _sparsePattern_index, _sizeof_sparsePattern_index);


  _sparsePattern_maxColors = _system->getA_sparsePattern_maxColors();
  */
}

const int Peer::reportErrorMessage(ostream& messageStream)
{
  /*
  if (_solverStatus == ISolver::SOLVERERROR)
  {
    if (_idid == -1)
      messageStream << "Invalid system dimension." << std::endl;
    if (_idid == -2)
      messageStream << "Method not implemented." << std::endl;
    if (_idid == -3)
      messageStream << "No valid system/settings available." << std::endl;
    if (_idid == -11)
      messageStream << "Step size too small." << std::endl;
  }

  else if (_solverStatus == ISolver::USER_STOP)
  {
    messageStream << "Simulation terminated by user at t: " << _tCurrent << std::endl;
  }

  return _idid;
  */
  return 0;
}

void Peer::writeSimulationInfo()
{
  /*
#ifdef USE_BOOST_LOG
  src::logger lg;

  // Now, let's try logging with severity
  src::severity_logger<cvodeseverity_level> slg;

  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nfQSe, netfQS;

  int qlast, qcur;
  realtype h0u, hlast, hcur, tcur;

  int flag;

  flag = CVodeGetIntegratorStats(_cvodeMem, &nst, &nfe, &nsetups, &netf, &qlast, &qcur, &h0u, &hlast, &hcur, &tcur);

  flag = CVodeGetNonlinSolvStats(_cvodeMem, &nni, &ncfn);

  BOOST_LOG_SEV(slg, cvode_normal)<< " Number steps: " << nst;
  BOOST_LOG_SEV(slg, cvode_normal)<< " Function evaluations " << "f: " << nfe;
  BOOST_LOG_SEV(slg, cvode_normal)<< " Error test failures " << "netf: " << netfS;
  BOOST_LOG_SEV(slg, cvode_normal)<< " Linear solver setups " << "nsetups: " << nsetups;
  BOOST_LOG_SEV(slg, cvode_normal)<< " Nonlinear iterations " << "nni: " << nni;
  BOOST_LOG_SEV(slg, cvode_normal)<< " Convergence failures " << "ncfn: " << ncfn;

#endif
  //// Solver
  //outputStream  << "\nSolver: " << getName()
  //  << "\nVerfahren: ";
  //if(_cvodesettings->iMethod == EulerSettings::EULERFORWARD)
  //  outputStream << " Expliziter Peer\n\n";
  //else if(_cvodesettings->iMethod == EulerSettings::EULERBACKWARD)
  //  outputStream << " Impliziter Peer\n\n";

  //// System
  //outputStream
  //  << "Dimension  des Systems (ODE):             " << (int)_dimSys << "\n";
  //// Status, Anzahl Schritte, Nullstellenzeugs
  //SolverDefaultImplementation::writeSimulationInfo(outputStream);

  //// Nullstellensuche
  //if (_cvodesettings->iZeroSearchMethod == SolverSettings::NO_ZERO_SEARCH)
  //{
  //  outputStream << "Nullstellensuche:                         Keine\n\n" << endl;
  //}
  //else
  //{
  //  if (_cvodesettings->iZeroSearchMethod == SolverSettings::BISECTION)
  //  {
  //  outputStream << "Nullstellensuche:                         Bisektion\n" << endl;
  //  }
  //  else
  //  {
  //  outputStream << "Nullstellensuche:                         Lineare Interpolation\n" << endl;
  //  }

  //}

  //// Schritteweite
  //outputStream
  //  << "ausgegebene Schritte:                     " << _outStps << "\n"
  //  << "Anfangsschrittweite:                      " << _cvodesettings->dH_init << "\n"
  //  << "Ausgabeschrittweite:                      " << dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings()->gethOutput() << "\n"
  //  << "Obere Grenze fÃ¼r Schrittweite:            " << _hUpLim << "\n\n";
  //// Status
  //outputStream
  //  << "Solver-Status:                            " << _idid << "\n\n";
  */
}

int Peer::check_flag(void *flagvalue, const char *funcname, int opt)
{
  /*
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated

  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return (1);
  }

  //Check if flag < 0

  else if (opt == 1)
  {
    errflag = (int *) flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return (1);
    }
  }

  // Check if function returned NULL pointer - no memory allocated

  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return (1);
  }

  return (0);
  */
  return 0;
}

void Peer::setcycletime(double cycletime){}
