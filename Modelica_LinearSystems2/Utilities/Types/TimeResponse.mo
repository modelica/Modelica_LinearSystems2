within Modelica_LinearSystems2.Utilities.Types;
type TimeResponse = enumeration(
    Impulse "Impulse response",
    Step "Step response",
    Ramp "Ramp response",
    Initial "Initial condition response") "Enumeration of time response type"
    annotation (Evaluate=true,
    Icon(graphics={Ellipse(
        extent={{-100,100},{100,-100}},
        lineColor={255,0,128},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid), Text(
        extent={{-94,94},{94,-94}},
        lineColor={255,0,128},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid,
        textString="e")}));
