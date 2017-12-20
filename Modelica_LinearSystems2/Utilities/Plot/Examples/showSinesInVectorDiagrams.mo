within Modelica_LinearSystems2.Utilities.Plot.Examples;
function showSinesInVectorDiagrams
  "Plot sine functions in several diagrams in vector layout"
  extends Modelica.Icons.Function;

  input Modelica.SIunits.Frequency freqHz1 = 2 "Frequency of sine wave 1";
  input Modelica.SIunits.Frequency freqHz2 = 3 "Frequency of sine wave 2";
  input Modelica.SIunits.Damping damping1 = 0.8
    "Damping coefficient of sine wave 12";
  input Modelica.SIunits.Damping damping2 = 0.1
    "Damping coefficient of sine wave 2";

protected
  Integer nX = 500;
  Integer nSymbol = 50;
  Integer nPeriod = 5;
  Real x1[nX];
  Real y1[nX];
  Real x2[nX];
  Real y2[nX];
  Real x3[nSymbol];
  Real y3[nSymbol];
algorithm
  (x1,y1) :=Utilities.dampedSine(freqHz1, damping1, nPeriod, nX);
  (x2,y2) :=Utilities.dampedSine(freqHz2, damping2, nPeriod, nX);
  (x3,y3) :=Utilities.dampedSine(freqHz1, 2*damping1, nPeriod, nSymbol);

  Modelica_LinearSystems2.Utilities.Plot.diagramVector({
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x1,
          y=y1,
          legend="torque1"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x2,
          y=y2,
          legend="torque2",
          autoLine=false,
          lineColor={255,0,0})},
      heading="Bearing friction torques",
      yLabel="[N.m]"),
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x1,
          y=1.1*y1,
          legend="torque3"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x2,
          y=1.2*y2,
          legend="torque4",
          autoLine=false,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Dot),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x3,
          y=y3,
          legend="torque5",
          autoLine=false,
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Circle)},
      xLabel="w [rad/s]",
      yLabel="[N.m]")});

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This function plots the following diagram (with default input arguments):
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showSinesInVectorDiagrams.png\">
</blockquote>
</html>"));
end showSinesInVectorDiagrams;
