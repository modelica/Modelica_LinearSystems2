within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode1
  "Construct two transfer functions and plot the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  output Boolean ok;

protected
  TransferFunction s = TransferFunction.s();

  TransferFunction tf1= 1/(10*s+1)^3;
  TransferFunction tf2=(s + 2)/(2*s^2 + 3*s +4);
algorithm
  TransferFunction.Plot.bode(
    tf=tf1);
  TransferFunction.Plot.bode(
    tf=tf2);
  ok := true;

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Construct transfer functions
</p>
<blockquote>
tf<sub>1</sub> = 1 / (10*<var>s</var> + 1)<sup>3</sup>
</blockquote>
<p>
and
</p>
<blockquote>
tf<sub>2</sub> = (<var>s</var> + 2) / (2*<var>s</var><sup>2</sup> + 3*<var>s</var> + 4)
</blockquote>
<p>
and plot their Bode diagram with automatic determination of the frequency range to plot.
</p>
</html>"));
end plotBode1;
