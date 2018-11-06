within Modelica_LinearSystems2.Utilities.Plot.Examples;
function showPointSymbols "Show the available point symbols"
  extends Modelica.Icons.Function;

protected
  Real y1[:]=0:0.2:1.0;
  Real x1[:]=0.01*ones(size(y1, 1));

  Real x2[:]=x1 + 0.1*ones(size(y1, 1));
  Real x3[:]=x2 + 0.1*ones(size(y1, 1));
  Real x4[:]=x3 + 0.1*ones(size(y1, 1));
  Real x5[:]=x4 + 0.1*ones(size(y1, 1));
  Real x6[:]=x5 + 0.1*ones(size(y1, 1));
  Real x7[:]=x6 + 0.1*ones(size(y1, 1));
  Real x8[:]=x7 + 0.1*ones(size(y1, 1));
algorithm
  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x1,
          y=y1,
          autoLine=false,
          lineColor={255,128,0},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Cross,
          legend="Cross"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x2,
          y=y1,
          autoLine=false,
          lineColor={0,128,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Circle,
          legend="Circle"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x3,
          y=y1,
          autoLine=false,
          lineColor={0,0,0},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Square,
          legend="Square"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x4,
          y=y1,
          autoLine=false,
          lineColor={255,0,0},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.FilledSquare,
          legend="FilledSquare"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x5,
          y=y1,
          autoLine=false,
          lineColor={0,0,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.TriangleDown,
          legend="TriangleDown"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x6,
          y=y1,
          autoLine=false,
          lineColor={255,0,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.TriangleUp,
          legend="TriangleUp"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x7,
          y=y1,
          autoLine=false,
          lineColor={128,255,128},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Diamond,
          legend="Diamond"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x8,
          y=y1,
          autoLine=false,
          lineColor={128,128,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Dot,
          legend="Dot")},
      heading="PointSymbols / lineSymbols"));

  annotation (
    Documentation(info="<html>
<p>
This function demonstrates the supported point symbols defined via enumeration
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol\">Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol</a>:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showPointSymbols.png\"></p>
</html>"));
end showPointSymbols;
