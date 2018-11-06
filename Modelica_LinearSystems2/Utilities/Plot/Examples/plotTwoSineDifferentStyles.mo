within Modelica_LinearSystems2.Utilities.Plot.Examples;
function plotTwoSineDifferentStyles
  "Plot two sine functions in one diagram with different styles"
  extends Modelica.Icons.Function;

  input Modelica.SIunits.Frequency freqHz1 = 2 "Frequency of sine wave 1";
  input Modelica.SIunits.Frequency freqHz2 = 3 "Frequency of sine wave 2";
  input Modelica.SIunits.Damping damping1 = 0.8
    "Damping coefficient of sine wave 12";
  input Modelica.SIunits.Damping damping2 = 0.1
    "Damping coefficient of sine wave 2";

protected
  Integer nX = 500;
  Integer nSymbol = 20;
  Integer nPeriod = 5;
  Real x1[nX];
  Real y1[nX];
  Real x2[nX];
  Real y2[nX];
algorithm
  (x1,y1) :=Utilities.dampedSine(freqHz1, damping1, nPeriod, nX);
  (x2,y2) :=Utilities.dampedSine(freqHz2, damping2, nPeriod, nX);

  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
      curve={
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x1,
          y=y1,
          legend="torque1",
          autoLine=false,
          lineColor={0,127,0},
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Circle),
        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=x2,
          y=y2,
          legend="torque2",
          autoLine=false,
          lineColor={255,0,0},
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Dot,
          lineThickness=0.75)},
      heading="Bearing friction torques",
      xLabel="w [rad/s]",
      yLabel="[N.m]"));

  annotation (
    Documentation(info="<html>
<p>
This function plots the following diagram (with default input arguments):
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotTwoSineDifferentStyles.png\">
</blockquote>
</html>"));
end plotTwoSineDifferentStyles;
