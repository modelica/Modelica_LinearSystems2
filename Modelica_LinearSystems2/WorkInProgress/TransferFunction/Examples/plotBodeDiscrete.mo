within Modelica_LinearSystems2.WorkInProgress.TransferFunction.Examples;
function plotBodeDiscrete
  "Plot the Bode diagram of the continuous and the discrete transfer functions with automatic determination of the frequency range to plot"

  import Modelica.Utilities.Streams.print;
  import Complex;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

  output Boolean ok;
protected
  Modelica.SIunits.Time Ts = 0.1 "Sample time";
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact "Discretization method";
  TransferFunction tf=TransferFunction({1}, {1,0,1});
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=
    Modelica_LinearSystems2.DiscreteTransferFunction(tf,Ts,method);

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.bode(tf);
  Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction.Plot.bode(dtf);
  ok := true;

end plotBodeDiscrete;
