within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBodeFilter3 "Show high pass filters of all filter types"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import AF = Modelica_LinearSystems2.Utilities.Types.AnalogFilter;
  import FT = Modelica_LinearSystems2.Utilities.Types.FilterType;

  input Integer order(min=1) = 4 "Order of filter";
  input Modelica.Units.SI.Frequency f_cut=1 "Cut-off frequency";
  input Real A_ripple(unit="dB") = 3
    "Pass band ripple for Chebyshev filter (otherwise not used)";
  input Modelica.Units.SI.Frequency f_min(min=0) = 0.1 "Minimum frequency value";
  input Modelica.Units.SI.Frequency f_max(min=0) = 10 "Maximum frequency value";
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
  ZerosAndPoles.Plot.bode(
    tf1,
    autoRange=false,
    f_min=f_min,
    f_max=f_max,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(
      heading="Bode plot: " + String(tf1) + ", critical damping filter",
      heightRatio=0.5));
  ZerosAndPoles.Plot.bode(
    tf2,
    autoRange=false,
    f_min=f_min,
    f_max=f_max,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(
      heading="Bode plot: " + String(tf2) + ", Bessel filter",
      heightRatio=0.5));
  ZerosAndPoles.Plot.bode(
    tf3,
    autoRange=false,
    f_min=f_min,
    f_max=f_max,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(
      heading="Bode plot: " + String(tf3) + ", Butterworth filter",
      heightRatio=0.5));
  ZerosAndPoles.Plot.bode(
    tf4,
    autoRange=false,
    f_min=f_min,
    f_max=f_max,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(
      heading="Bode plot: " + String(tf4) + ", Chebyshev filter",
      heightRatio=0.5));

  ok := true;

  annotation (__Dymola_interactive=true);
end plotBodeFilter3;
