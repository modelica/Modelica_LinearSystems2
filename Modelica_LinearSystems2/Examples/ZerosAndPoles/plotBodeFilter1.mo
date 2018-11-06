within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBodeFilter1 "Example for computation of filter and plot frequency response of filter"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Utilities.Types;

  input Utilities.Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
  input Integer order=2;
  input Modelica.SIunits.Frequency f_cut=10;
  output Boolean ok;
protected
  ZerosAndPoles tf_filter=ZerosAndPoles.Design.filter(
    order=order,
    f_cut=f_cut,
    analogFilter=analogFilter);

algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(tf_filter);
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows how to construct a ZerosAndPoles transfer function description of a filter
and to plot frequency response of such a filter.
</p>
</html>
"));
end plotBodeFilter1;
