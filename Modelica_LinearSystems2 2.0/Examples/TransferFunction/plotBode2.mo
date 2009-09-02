within Modelica_LinearSystems2.Examples.TransferFunction;
function plotBode2
  "Construct 2 transfer functions and plot the Bode diagram with automatic determination of the frequency range to plot"

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.TransferFunction;

  output Boolean ok;
  annotation (interactive=true);
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

end plotBode2;
