#pragma once
/** @addtogroup simcorefactoriesPolicies
 *
 *  @{
 */

boost::shared_ptr<ISimData> createSimDataFunction();
boost::shared_ptr<ISimVars> createSimVarsFunction(size_t dim_real,size_t dim_int,size_t dim_bool,size_t dim_string,size_t dim_pre_vars,size_t dim_z,size_t z_i);
boost::shared_ptr<IMixedSystem> createSystemFunction(IGlobalSettings* globalSettings,boost::shared_ptr<IAlgLoopSolverFactory> algloopsolverfactory,boost::shared_ptr<ISimData> simData,boost::shared_ptr<ISimVars> simVars);
boost::shared_ptr<IAlgLoopSolverFactory> createStaticAlgLoopSolverFactory(IGlobalSettings* globalSettings,PATH library_path,PATH modelicasystem_path);
/*
Policy class to create a OMC-,  Modelica- system or AlgLoopSolver
*/
template <class CreationPolicy>
struct StaticSystemOMCFactory: public  ObjectFactory<CreationPolicy>
{
public:
  StaticSystemOMCFactory(PATH library_path, PATH modelicasystem_path, PATH config_path)
    :ObjectFactory<CreationPolicy>(library_path, modelicasystem_path, config_path)
  {
      _use_modelica_compiler = false;
  }

  virtual ~StaticSystemOMCFactory()
  {
  }

  virtual boost::shared_ptr<ISimData> createSimData()
  {
    return createSimDataFunction();
  }

  virtual boost::shared_ptr<ISimVars> createSimVars(size_t dim_real,size_t dim_int,size_t dim_bool,size_t dim_string,size_t dim_pre_vars,size_t dim_z,size_t z_i)
  {
    return createSimVarsFunction(dim_real, dim_int, dim_bool, dim_string, dim_pre_vars, dim_z, z_i);
  }

  virtual boost::shared_ptr<IAlgLoopSolverFactory> createAlgLoopSolverFactory(IGlobalSettings* globalSettings)
  {
    return createStaticAlgLoopSolverFactory(globalSettings,ObjectFactory<CreationPolicy>::_library_path,ObjectFactory<CreationPolicy>::_modelicasystem_path);
  }

  boost::shared_ptr<IMixedSystem> createSystem(string modelLib,string modelKey,IGlobalSettings* globalSettings,boost::shared_ptr<IAlgLoopSolverFactory> algloopsolverfactory,boost::shared_ptr<ISimData> simData,boost::shared_ptr<ISimVars> simVars)
  {
     return createSystemFunction(globalSettings,algloopsolverfactory,simData,simVars);
  }
   boost::shared_ptr<IMixedSystem>  createModelicaSystem(PATH modelica_path, string modelKey, IGlobalSettings* globalSettings, boost::shared_ptr<IAlgLoopSolverFactory> algloopsolverfactory,boost::shared_ptr<ISimData> simData,boost::shared_ptr<ISimVars> simVars)
  {
    throw ModelicaSimulationError(MODEL_FACTORY,"Modelica is not supported");
  }
   bool _use_modelica_compiler;
};
/** @} */ // end of simcorefactoriesPolicies
