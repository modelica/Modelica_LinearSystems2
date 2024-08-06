within Modelica_LinearSystems2.Utilities.Plot.Examples;
function plotTwoSine "Plot two sine functions in one diagram"
  extends Modelica.Icons.Function;

  input Modelica.Units.SI.Frequency freqHz1=2 "Frequency of sine wave 1";
  input Modelica.Units.SI.Frequency freqHz2=3 "Frequency of sine wave 2";
  input Modelica.Units.SI.Damping damping1=0.8
    "Damping coefficient of sine wave 12";
  input Modelica.Units.SI.Damping damping2=0.1
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
    curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=x1,
      y=y1,
      legend="torque1"),Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=x2,
      y=y2,
      legend="torque2")},
    heading="Bearing friction torques",
    xLabel="w [rad/s]",
    yLabel="[N.m]"));
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots the following diagram (with default input arguments):
</p>

<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotTwoSine.png\">
</div>
</html>"));
end plotTwoSine;
