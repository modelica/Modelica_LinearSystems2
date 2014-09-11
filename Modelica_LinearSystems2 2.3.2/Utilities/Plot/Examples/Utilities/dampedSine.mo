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

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(x,y) = Plot.Examples.Utilities.<b>dampedSine</b>(freqHz, damping, nPeriod=5, nPoints=500);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes nPoints points x[i], y[i] of a sine-curve
</p>

<pre>
   y = exp(-x*damping)*sin(2*pi*freqHz*x)
</pre>

<p>
where the abszissa values x[i] are in the range 0 .. nPeriod*T, where T = 1/freqHz is the period
of the sine wave.
</p>

<h4>Example</h4>
<p>
With the default options and freqHz=2 Hz and damping=0.8, the following curve is generated:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/dampedSine.png\"/></p>
</html>"));
end dampedSine;
