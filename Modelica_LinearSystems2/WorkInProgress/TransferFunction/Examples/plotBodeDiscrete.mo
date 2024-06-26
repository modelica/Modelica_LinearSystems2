within Modelica_LinearSystems2.WorkInProgress.TransferFunction.Examples;
function plotBodeDiscrete
  "Plot the Bode diagram of the continuous and the discrete transfer functions with automatic determination of the frequency range to plot"

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

  output Boolean ok;
protected
  Modelica.Units.SI.Time Ts=0.1 "Sample time";
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact "Discretization method";
  TransferFunction tf=TransferFunction({1}, {1,0,1});
  DiscreteTransferFunction dtf=DiscreteTransferFunction(tf,Ts,method);

algorithm
  TransferFunction.Plot.bode(tf);
  DiscreteTransferFunction.Plot.bode(dtf);
  ok := true;

  annotation (__Dymola_interactive=true);
end plotBodeDiscrete;
