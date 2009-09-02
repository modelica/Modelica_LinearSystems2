within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode1
  "Construct a transfer function and plot the Bode diagram with automatic determination of the frequency range to plot"

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.TransferFunction;

  output Boolean ok;

  annotation (interactive=true);

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

end plotBode1;
