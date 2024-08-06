within Modelica_LinearSystems2.Controllers.Internal;
function convertToInit "Convert InitWithGlobalDefault to Init"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Controllers.Types;

  input Types.InitWithGlobalDefault globalInit;
  input Types.Init sampleInit;
  output Types.Init Init;
algorithm
  Init := if globalInit == Types.InitWithGlobalDefault.UseSampleClockOption then
                sampleInit
            else if globalInit == Types.InitWithGlobalDefault.NoInit then
                Types.Init.NoInit
            else if globalInit == Types.InitWithGlobalDefault.SteadyState then
                Types.Init.SteadyState
            else if globalInit == Types.InitWithGlobalDefault.InitialState then
                Types.Init.InitialState
            else Types.Init.InitialOutput;
end convertToInit;
