within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode1
  "Construct a transfer function and plot the Bode diagram with automatic determination of the frequency range to plot"
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

  annotation (__Dymola_interactive=true);
end plotBode1;
