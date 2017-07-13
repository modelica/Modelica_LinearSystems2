within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode2
  "Construct 2 transfer functions and plot the Bode diagram with automatic determination of the frequency range to plot"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Complex;

  output Boolean ok;
protected
  TransferFunction tf1=TransferFunction({1,2}, {2,3,4});
  TransferFunction tf2=TransferFunction({1,2}, {4,1,4});

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.bode(
    tf1,
    autoRange=false,
    f_min=0.01,
    f_max=30);
  Modelica_LinearSystems2.TransferFunction.Plot.bode(
    tf2,
    autoRange=false,
    f_min=0.01,
    f_max=30);
  ok := true;

  annotation (__Dymola_interactive=true);
end plotBode2;
