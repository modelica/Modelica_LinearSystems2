within Modelica_LinearSystems2.Utilities.Plot.Examples;
function showLegendStyles
  "Show several vector-diagram plots that demonstrate the various legend options"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Utilities.Plot.Records;
  import Modelica_LinearSystems2.Utilities.Plot.Types;

protected
  Real x[2]={0,1};
  Real y1[2]={0,0.5};
  Real y2[2]={0.2,0.7};
  Real y3[2]={0.4,0.9};
  Records.Curve curves[3]={
    Records.Curve(x=x, y=y1, legend="curve 1"),
    Records.Curve(x=x, y=y2, legend="curve 2"),
    Records.Curve(x=x, y=y3, legend="curve 3")};

algorithm
  Modelica_LinearSystems2.Utilities.Plot.diagramVector({
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.Above",
        legendLocation=Types.LegendLocation.Above),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.Above, legendHorizontal=false",
        legendHorizontal=false,
        legendLocation=Types.LegendLocation.Above),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.Right",
        legendLocation=Types.LegendLocation.Right),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.Below",
        legendLocation=Types.LegendLocation.Below)},
    Records.Device(
      xTopLeft=0,
      yTopLeft=0,
      diagramWidth=100));

  Modelica_LinearSystems2.Utilities.Plot.diagramVector({
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.TopLeft",
        legendLocation=Types.LegendLocation.TopLeft),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.TopRight, legendHorizontal=false",
        legendHorizontal=false,
        legendLocation=Types.LegendLocation.TopRight),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.BottomLeft",
        legendLocation=Types.LegendLocation.BottomLeft),
      Records.Diagram(
        curve=curves,
        heading="LegendLocation.BottomRight",
        legendLocation=Types.LegendLocation.BottomRight)},
    Records.Device(
      xTopLeft=105,
      yTopLeft=0,
      diagramWidth=100));

  annotation(
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This function plots the following diagram (with default input arguments):
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showLegendStyles.png\">
</blockquote>
</html>"));
end showLegendStyles;
