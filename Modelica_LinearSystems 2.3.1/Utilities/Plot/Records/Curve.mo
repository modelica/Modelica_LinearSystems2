within Modelica_LinearSystems2.Utilities.Plot.Records;
record Curve "Properties of a curve (displayed in a diagram)"
  extends Modelica.Icons.Record;

  Real x[:] "Values x of curve" annotation(Dialog);
  Real y[:] "Values y of curve" annotation(Dialog);
  String legend="" "Legend text of curve" annotation(Dialog);

  Boolean autoLine = true "True, if automatic line properties of curve"
    annotation(  choices(checkBox=true));

  Integer lineColor[3]={0,0,255} "Color of curve as rgb values"
    annotation(Dialog(enable=not autoLine,group="If autoLine = false (otherwise ignored)",colorSelector=true));

  Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern linePattern=
    Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid
    "Line pattern of curve" annotation(Dialog(enable=not autoLine,group="If autoLine = false (otherwise ignored)"));
  Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol lineSymbol=
    Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None
    "Symbol for points on curve" annotation(Dialog(enable=not autoLine,group="If autoLine = false (otherwise ignored)"));
  Real lineThickness=0.25 "Line thickness of curve"
    annotation(Dialog(group="If autoLine = false (otherwise ignored)"));

/*
  Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineSymbolSize=3
    "Symbol size" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
*/

  annotation (Documentation(info="<html>
<p>
With this record the properties of one curve in a diagram are defined
(such as data points, color, line pattern).
</p>
</html>"));
end Curve;
