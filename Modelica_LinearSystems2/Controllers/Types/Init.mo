within Modelica_LinearSystems2.Controllers.Types;
type Init = enumeration(
    NoInit
      "No initialization (start values are used as guess values with fixed=false)",
    SteadyState "Steady state initialization (derivatives of states are zero)",
    InitialState "Initialization with initial states",
    InitialOutput
      "Initialization with initial outputs (and steady state of the states if possibles)")
  "Enumeration defining initialization of a block"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));
