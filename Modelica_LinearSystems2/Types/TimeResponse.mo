within Modelica_LinearSystems2.Types;
type TimeResponse = enumeration(
    Impulse "Impulse response",
    Step "Step response",
    Ramp "Ramp response",
    Initial "Initial condition response") "Enumeration of time response type"
    annotation (
      Evaluate=true,
      obsolete = "Obsolete enumeration - use Modelica_LinearSystems2.Utilities.Types.TimeResponse instead");
