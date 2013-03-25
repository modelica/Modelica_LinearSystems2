within Modelica_LinearSystems2.Utilities.Plot.Records;
record RootLocusDiagram "Properties of a root locus diagram"
  extends Modelica.Icons.Record;

  String heading=""
    "Heading displayed above diagram (if empty, default heading)" annotation(Dialog);
  String ReName="Re"
    "Name of the real part of the eigen values (shown in tool tip)" annotation(Dialog);
  String ImName="Im"
    "Name of the imaginary part of the eigen values (shown in tool tip)" annotation(Dialog);
  Real heightRatio = 0.8 "Height of diagram = heightRatio*diagramWidth" annotation(Dialog);
  Boolean grid=true "True, if grid is shown" annotation(Dialog,  choices(__Dymola_checkBox=true));
  Boolean labelWithParam=false
    "True, if values of parameter shall be shown along the curves"
    annotation(Dialog,  choices(__Dymola_checkBox=true));

  /* group "Axes" (Axes properties) */
  String xLabel="Real part of eigenvalues"
    "String displayed at horizontal axis" annotation(Dialog(group="Axes"));
  String yLabel="Imaginary part of eigenvalues"
    "String displayed at vertical axis" annotation(Dialog(group="Axes"));
  Boolean logX = false "True, if logarithmic scale of x-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
  Boolean logY = false "True, if logarithmic scale of y-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
  Boolean uniformScaling = false
    "True, if same vertical and horizontal axis increment"
      annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));

   Integer lineColor[3]={0,0,255} "Color of curve as rgb values"
     annotation(Dialog(group="Curve properties",__Dymola_colorSelector, __Dymola_treeView=false));

   Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern linePattern=
      Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None
    "Line pattern of curve" annotation(Dialog(group="Curve properties"));
   Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol lineSymbol=
      Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Dot
    "Symbol for points on curve" annotation(Dialog(group="Curve properties"));
   Real lineThickness=0.25 "Line thickness of curve"
                              annotation(Dialog(group="Curve properties"));

  annotation (Documentation(info="<html>
<p>
With this record the properties of a root locus diagram can be defined, as needed by
function
<a href=\"Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel\">Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel</a>.
</p>
</html>"));
end RootLocusDiagram;
