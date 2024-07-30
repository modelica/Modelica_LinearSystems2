within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode3
  "Construct a transfer function from numerator and denominator array and plot the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  output Boolean ok;
protected
  TransferFunction tf=TransferFunction({1}, {1,0,1});

algorithm
  TransferFunction.Plot.bode(tf);
  ok := true;

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Construct a&nbsp;transfer function
</p>
<blockquote>
tf = 1 / (<var>s</var><sup>2</sup> + 1)
</blockquote>
<p>
from numerator and denominator array and plot the Bode diagram with automatic determination
of the frequency range to plot.
</p>
</html>"));
end plotBode3;
