within Modelica_LinearSystems2.Utilities.Types;
type Window = enumeration(
    Rectangle "Rectangle",
    Bartlett "Bartlett",
    Hann "Hann",
    Hamming "Hamming",
    Blackman "Blackman",
    Kaiser "Kaiser") "Enumeration of window type"
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
