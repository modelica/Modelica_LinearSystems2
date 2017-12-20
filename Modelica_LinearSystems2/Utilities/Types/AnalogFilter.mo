within Modelica_LinearSystems2.Utilities.Types;
type AnalogFilter = enumeration(
    CriticalDamping "Filter with critical damping",
    Bessel "Bessel filter",
    Butterworth "Butterworth filter",
    Chebyshev "Chebyshev filter")
  "Enumeration defining the method of filtering"
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
