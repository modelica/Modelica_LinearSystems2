within Modelica_LinearSystems2.Utilities.Plot.Records;
record CurveProperties "Properties of a curve"
  extends Modelica.Icons.Record;

   Integer lineColor[3]={0,0,255} "Color of curve as rgb values"
     annotation(Dialog(group="Curve properties",__Dymola_colorSelector, __Dymola_treeView=false));

   Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern linePattern=
      Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid
    "Line pattern of curve" annotation(Dialog(group="Curve properties"));
   Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol lineSymbol=
      Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None
    "Symbol for points on curve" annotation(Dialog(group="Curve properties"));
   Real lineThickness=0.25 "Line thickness of curve"
                              annotation(Dialog(group="Curve properties"));

/*
   Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineSymbolSize=3
    "Symbol size" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
*/

end CurveProperties;
