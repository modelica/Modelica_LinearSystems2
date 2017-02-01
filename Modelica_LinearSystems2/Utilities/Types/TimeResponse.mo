within Modelica_LinearSystems2.Utilities.Types;
type TimeResponse = enumeration(
    Impulse "Impulse response",
    Step "Step response",
    Ramp "Ramp response",
    Initial "Initial condition response") "Enumeration of time response type"
    annotation (Evaluate=true, Documentation(info="<html>
</html>"));
