within Modelica_LinearSystems2.Internal;
function methodString
  "Return a sting with the name of the discretization method"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Utilities.Types.Method;

  input Method method;
  output String s;
algorithm
s := if method==Method.ExplicitEuler then "ExplicitEuler" else
          if method==Method.ImplicitEuler then "ImplicitEuler" else
          if method==Method.Trapezoidal then "Trapezoidal" else
          if method==Method.ImpulseExact then "ImpulseExact" else
          if method==Method.StepExact then "StepExact" else
          "RampExact";
end methodString;
