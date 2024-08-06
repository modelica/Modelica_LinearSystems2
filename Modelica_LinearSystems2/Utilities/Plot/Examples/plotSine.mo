within Modelica_LinearSystems2.Utilities.Plot.Examples;
function plotSine "Plot a sine function in one diagram"
  extends Modelica.Icons.Function;

  input Modelica.Units.SI.Frequency freqHz=2 "Frequency of sine wave";
  input Modelica.Units.SI.Damping damping=0.8
    "Damping coefficient of sine wave";
protected
   Integer nX = 500;
   Integer nPeriod = 5;
   Real x[nX];
   Real y[nX];
algorithm
   (x,y) :=Utilities.dampedSine(freqHz, damping, nPeriod, nX);

  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
    curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
      x=x,
      y=y,
      legend="torque1")},
    heading="Bearing friction torque",
    xLabel="w [rad/s]",
    yLabel="[N.m]"));
    annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots the following diagram (with default input arguments):
</p>
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotSine.png\">
</div>
</html>"));
end plotSine;
