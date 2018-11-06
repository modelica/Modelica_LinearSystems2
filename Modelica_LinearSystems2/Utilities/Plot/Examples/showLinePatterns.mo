within Modelica_LinearSystems2.Utilities.Plot.Examples;
function showLinePatterns "Show the available line patterns"
  extends Modelica.Icons.Function;

protected
  Real x1[:]={0,1};
  Real y1[:]=x1;

  Real x2[:]=x1 + 0.1*ones(size(x1, 1));
  Real x3[:]=x2 + 0.1*ones(size(x1, 1));
  Real x4[:]=x3 + 0.1*ones(size(x1, 1));
  Real x5[:]=x4 + 0.1*ones(size(x1, 1));
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
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
          legend="Solid"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x2,
          y=y1,
          autoLine=false,
          lineColor={0,128,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Dash,
          legend="Dash"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x3,
          y=y1,
          autoLine=false,
          lineColor={0,0,0},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Dot,
          legend="Dot"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x4,
          y=y1,
          autoLine=false,
          lineColor={255,0,0},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.DashDot,
          legend="DashDot"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x5,
          y=y1,
          autoLine=false,
          lineColor={0,0,255},
          lineThickness=0.5,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.DashDotDot,
          legend="DashDotDot")},
      heading="Line patterns"));

  annotation (
    Documentation(info="<html>
<p>
This function demonstrates the supported line patters defined via enumeration
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern\">Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern</a>:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showLinePatterns.png\"></p>
</html>"));
end showLinePatterns;
