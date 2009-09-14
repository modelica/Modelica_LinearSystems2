within Modelica_LinearSystems2.Types;
type Method = enumeration(
    ExplicitEuler "discretize block using explicit Euler integration",
    ImplicitEuler "discretize block using implicit Euler integration",
    Trapezoidal
      "discretize block using trapezoidal integration (Tustins method)",
    StepExact "discretize block using the zero-order hold equivalent",
    RampExact "discretize block using first-order hold equivalent",
    ImpulseExact "discretize block using impulse invariant discretization")
  "Enumeration of discretization methods" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));
