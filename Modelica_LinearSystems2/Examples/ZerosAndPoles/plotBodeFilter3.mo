within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBodeFilter3 "Example showing high-pass filters of all filter types"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import AF = Modelica_LinearSystems2.Utilities.Types.AnalogFilter;
  import FT = Modelica_LinearSystems2.Utilities.Types.FilterType;

  input Integer order(min=1) = 4 "Order of filter";
  input Modelica.SIunits.Frequency f_cut=1 "Cut-off frequency";
  input Real A_ripple(unit="dB") = 3
    "Pass band ripple for Chebyshev filter (otherwise not used)";
  input Modelica.SIunits.Frequency f_min(min=0) = 0.1 "Minimum frequency value";
  input Modelica.SIunits.Frequency f_max(min=0) = 10 "Maximum frequency value";
  output Boolean ok;
protected
  ZerosAndPoles tf1=ZerosAndPoles.Design.filter(
    analogFilter=AF.CriticalDamping,
    order=order,
    f_cut=f_cut,
    A_ripple=A_ripple,
    filterType=FT.HighPass);
  ZerosAndPoles tf2=ZerosAndPoles.Design.filter(
    analogFilter=AF.Bessel,
    order=order,
    f_cut=f_cut,
    A_ripple=A_ripple,
    filterType=FT.HighPass);
  ZerosAndPoles tf3=ZerosAndPoles.Design.filter(
    analogFilter=AF.Butterworth,
    order=order,
    f_cut=f_cut,
    A_ripple=A_ripple,
    filterType=FT.HighPass);
  ZerosAndPoles tf4=ZerosAndPoles.Design.filter(
    analogFilter=AF.Chebyshev,
    order=order,
    f_cut=f_cut,
    A_ripple=A_ripple,
    filterType=FT.HighPass);
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    tf1,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    tf2,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    tf3,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    tf4,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);

  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example constructs ZerosAndPoles transfer function descriptions of various <b>high-pass filters</b>
and plots correspondent frequency responses.
</p>
</html>
"));
end plotBodeFilter3;
