within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBodeFilter1 "Compute filter and plot frequency response of filter"
  import Modelica;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Types;

  input Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping
    "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
  input Integer order=2;
  input Modelica.SIunits.Frequency f_cut=10;
  output Boolean ok;
  annotation (interactive=true);
  annotation (interactive=true);
protected
  ZerosAndPoles tf_filter=ZerosAndPoles.Design.filter(
      order=order,
      f_cut=f_cut,
      analogFilter=analogFilter);

algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
                                  tf_filter);
  ok := true;
equation

end plotBodeFilter1;
