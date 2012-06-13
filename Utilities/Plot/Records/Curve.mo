within Modelica_LinearSystems2.Utilities.Plot.Records;
record Curve "Properties of a curve (displayed in a diagram)"
  extends Modelica.Icons.Record;

   Real x[:] "x-values of curve" annotation(Dialog);
   Real y[:] "y-values of curve" annotation(Dialog);
   String legend="" "Legend text of curve" annotation(Dialog);

   Boolean autoLine = true "True, if automatic line properties of curve"
     annotation(Dialog,  choices(__Dymola_checkBox=true));

   Integer lineColor[3]={0,0,255} "Color of curve as rgb values"
     annotation(Dialog(group="If autoLine = false (otherwise ignored)",__Dymola_colorSelector, __Dymola_treeView=false));

   Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern linePattern=
      Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid
    "Line pattern of curve" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
   Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol lineSymbol=
      Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None
    "Symbol for points on curve" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));

/*
   Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineThickness=0.25 
    "Line thickness of curve" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
   Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineSymbolSize=3 
    "Symbol size" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
*/

end Curve;
