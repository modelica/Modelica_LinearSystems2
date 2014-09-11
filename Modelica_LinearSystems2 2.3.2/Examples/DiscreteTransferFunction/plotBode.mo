within Modelica_LinearSystems2.Examples.DiscreteTransferFunction;
function plotBode "Bode plot of a DiscreteTransferFunction objekt"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

protected
  Real Ts=0.1;
  DiscreteTransferFunction z=DiscreteTransferFunction.z(Ts);
  DiscreteTransferFunction dtf=(0.000160362*z^5 + 0.000197116*z^4 - 0.0011001*z^3 + 0.00091796*z^2 - 7.25742e-005*z - 9.97232e-005)/(z^6 - 5.45248*z^5 + 12.4075*z^4 - 15.0825*z^3 + 10.3296*z^2 - 3.77899*z + 0.57695);

  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
algorithm
  dtf.Ts := Ts;
  dtf.method := method;
  DiscreteTransferFunction.Plot.bode(dtf);

annotation (__Dymola_interactive=true);
end plotBode;
