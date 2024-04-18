within Modelica_LinearSystems2.Controllers.Types;
type InitWithGlobalDefault = enumeration(
    NoInit
      "No initialization (start values are used as guess values with fixed=false)",
    SteadyState "Steady state initialization (derivatives of states are zero)",
    InitialState "Initialization with initial states",
    InitialOutput
      "Initialization with initial outputs (and steady state of the states if possibles)",
    UseSampleClockOption "Use initialization defined in sampleClock component")
  "Enumeration defining initialization of a block"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));
