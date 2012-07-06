within Modelica_LinearSystems2.Utilities.Plot.Records;
record RootLocusDiagram "Properties of a root locus diagram"
  extends Modelica.Icons.Record;

  String heading=""
    "Heading displayed above diagram (if empty, default heading)"                 annotation(Dialog);
  Boolean grid=true "True, if grid is shown" annotation(Dialog,  choices(__Dymola_checkBox=true));

/*
  Integer markerColorMin[3]={0,0,255}
    "Color of marker for minimum parameter value"
    annotation(Dialog(colorSelector=true));
  Integer markerColorMax[3]={255,0,0}
    "Color of marker for maximum parameter value"
    annotation(Dialog(colorSelector=true));
*/

  /* tab "Axes" (Axes properties) */
  String xLabel="Real part of eigenvalues"
    "String displayed at horizontal axis"                                        annotation(Dialog(tab="Axes"));
  String yLabel="Imaginary part of eigenvalues"
    "String displayed at vertical axis"                                             annotation(Dialog(tab="Axes"));
  Boolean logX = false "True, if logarithmic scale of x-axis" annotation(Dialog(tab="Axes"),choices(__Dymola_checkBox=true));
  Boolean logY = false "True, if logarithmic scale of y-axis" annotation(Dialog(tab="Axes"),choices(__Dymola_checkBox=true));
  Boolean uniformScaling = false
    "True, if same vertical and horizontal axis increment"
      annotation(Dialog(tab="Axes"),choices(__Dymola_checkBox=true));

  /* tab "Position and Size */
  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm xTopLeft=0
    "Horizontal position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog(tab="Position and Size"));
  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm yTopLeft=0
    "Vertical position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog(tab="Position and Size"));
  Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm diagramWidth=140
    "Width of diagram" annotation(Dialog(tab="Position and Size"));
  Real heightRatio = 0.7 "Height of diagram = heightRatio*diagramWidth" annotation(Dialog(tab="Position and Size"));
  Modelica_LinearSystems2.Utilities.Plot.Types.ImageResolution_dpi windowResolution=96
    "[dpi] Image resolution in window if applicable (e.g. unscaled window)"   annotation(Dialog(tab="Position and Size"));

end RootLocusDiagram;
