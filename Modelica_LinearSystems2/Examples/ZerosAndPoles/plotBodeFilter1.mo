within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBodeFilter1 "Compute filter and plot frequency response of filter"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Utilities.Types;

  input Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
  input Integer order=2;
  input Modelica.Units.SI.Frequency f_cut=10;
  output Boolean ok;
protected
  ZerosAndPoles tf_filter=ZerosAndPoles.Design.filter(
    order=order,
    f_cut=f_cut,
    analogFilter=analogFilter);

algorithm
  ZerosAndPoles.Plot.bode(tf_filter);
  ok := true;

  annotation (__Dymola_interactive=true);
end plotBodeFilter1;
