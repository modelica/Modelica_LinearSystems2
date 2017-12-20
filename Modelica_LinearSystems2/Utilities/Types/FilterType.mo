within Modelica_LinearSystems2.Utilities.Types;
type FilterType = enumeration(
    LowPass "Low pass filter",
    HighPass "High pass filter",
    BandPass "Band pass filter",
    BandStop "Band stop / notch filter")
  "Enumeration of analog filter types (high pass or low pass)"
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
