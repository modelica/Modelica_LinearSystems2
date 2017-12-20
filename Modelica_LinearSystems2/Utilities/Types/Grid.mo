within Modelica_LinearSystems2.Utilities.Types;
type Grid = enumeration(
    OneValue,
    Equidistant,
    Logarithmic) "Enumeration defining the grid type for a parameter variation"
  annotation (
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
