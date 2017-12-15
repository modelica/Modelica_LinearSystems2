within Modelica_LinearSystems2.Controllers.Types;
type MethodWithGlobalDefault = enumeration(
    ExplicitEuler "Discretization with explicit Euler integration",
    ImplicitEuler "Discretization with implicit Euler integration",
    Trapezoidal
      "Discretization with trapezoidal integration (Tustins method, recommended)",
    ImpulseExact "Exact discretization for impulse inputs",
    StepExact
      "Exact discretization for step inputs (zero-order hold equivalent)",
    RampExact
      "Exact discretization for ramp inputs (first-order hold equivalent)",
    UseSampleClockOption "Use method defined in sampleClock component")
  "Enumeration of discretization methods"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));
