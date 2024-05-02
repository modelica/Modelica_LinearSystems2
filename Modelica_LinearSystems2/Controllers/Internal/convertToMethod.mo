within Modelica_LinearSystems2.Controllers.Internal;
function convertToMethod "Convert MethodWithGlobalDefault to Method"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Controllers.Types;

  input Types.MethodWithGlobalDefault globalMethod;
  input Types.Method sampleMethod;
  output Types.Method method;
algorithm
  method := if globalMethod == Types.MethodWithGlobalDefault.UseSampleClockOption then
                sampleMethod
            else if globalMethod == Types.MethodWithGlobalDefault.ExplicitEuler then
                Types.Method.ExplicitEuler
            else if globalMethod == Types.MethodWithGlobalDefault.ImplicitEuler then
                Types.Method.ImplicitEuler
            else if globalMethod == Types.MethodWithGlobalDefault.Trapezoidal then
                Types.Method.Trapezoidal
            else if globalMethod == Types.MethodWithGlobalDefault.StepExact then
                Types.Method.StepExact
            else Types.Method.RampExact;
end convertToMethod;
