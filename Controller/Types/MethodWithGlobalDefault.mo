within Modelica_LinearSystems2.Controller.Types;
type MethodWithGlobalDefault = enumeration(
    ExplicitEuler "discretize block using explicit Euler integration",
    ImplicitEuler "discretize block using implicit Euler integration",
    Trapezoidal
      "discretize block using trapezoidal integration (Tustins method)",
    StepExact "discretize block using the zero-order hold equivalent",
    RampExact "discretize block using first-order hold equivalent",
    ImpulseExact "Impulse invariant discretization",
    UseSampleClockOption "use method defined in sampleClock component")
  "Enumeration of discretization methods" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));
