within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode1
  "Example for construction of two transfer functions and plot of the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Complex;

  output Boolean ok;

protected
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();

  TransferFunction tf1= 1/(10*s+1)^3;
  TransferFunction tf2=(s + 2)/(2*s^2 + 3*s +4);
algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.bode(
    tf=tf1);
  Modelica_LinearSystems2.TransferFunction.Plot.bode(
    tf=tf2);
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows how to construct a transfer function and to plot the Bode diagram
with automatic determination of the frequency range to plot.
</p>
</html>
"));
end plotBode1;
