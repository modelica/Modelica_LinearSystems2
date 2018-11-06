within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode3
  "Example for construction of transfer function and plot of the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Complex;

  output Boolean ok;
protected
  TransferFunction tf=TransferFunction({1}, {1,0,1});

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.bode(tf);
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows how to construct a transfer function and to plot the Bode diagram
with automatic determination of the frequency range to plot.
</p>
</html>
"));
end plotBode3;
