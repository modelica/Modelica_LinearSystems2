within Modelica_LinearSystems2.Controller.Types;
type MethodWithGlobalDefault = enumeration(
    ExplicitEuler "Discretization with explicit Euler integration",
    ImplicitEuler "Discretization with implicit Euler integration",
    Trapezoidal
      "Recommended: Discretization with trapezoidal integration (Tustins method)", 

    StepExact
      "Exact discretization for step inputs (zero-order hold equivalent)",
    RampExact
      "Exact discretization for ramp inputs (first-order hold equivalent)",
    ImpulseExact "Exact discretization for impulse inputs",
    UseSampleClockOption "Use method defined in sampleClock component")
  "Enumeration of discretization methods" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));

