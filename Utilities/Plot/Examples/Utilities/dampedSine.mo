within Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities;
function dampedSine "Return a damped sine chracteristic"
   input Modelica.SIunits.Frequency freqHz "Frequency of sine wave";
   input Modelica.SIunits.Damping damping "Damping coefficient of sine wave";
   input Integer nPeriod=5 "Number of periods to show";
   input Integer nPoints(min=2)=500 "Number of points";
   output Real x[nPoints];
   output Real y[nPoints];
protected
 Real xEnd = 1/freqHz*nPeriod;
algorithm
   x :=0:xEnd/(nPoints - 1):xEnd;
   for i in 1:size(x,1) loop
      y[i] :=Modelica.Math.exp(-x[i]*damping)*
             Modelica.Math.sin(2*Modelica.Constants.pi*freqHz*x[i]);
   end for;
end dampedSine;
