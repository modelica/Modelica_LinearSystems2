within Modelica_LinearSystems2.Types;
type Method = enumeration(
    ExplicitEuler "Discretization with explicit Euler integration",
    ImplicitEuler "Discretization with implicit Euler integration",
    Trapezoidal
      "Discretization with trapezoidal integration (Tustins method, recommended)", 

    ImpulseExact "Exact discretization for impulse inputs",
    StepExact
      "Exact discretization for step inputs (zero-order hold equivalent)",
    RampExact
      "Exact discretization for ramp inputs (first-order hold equivalent)")
  "Enumeration of discretization methods" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

