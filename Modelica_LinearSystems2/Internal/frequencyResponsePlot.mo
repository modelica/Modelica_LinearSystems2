within Modelica_LinearSystems2.Internal;
encapsulated function frequencyResponsePlot "Bode plot given f,A,phi values"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.Internal;
    import Modelica_LinearSystems2.Utilities.Plot;
    import SI = Modelica.SIunits;

  input Real f[:] "Frequency vector (either in Hz or rad/s)";
  input Real a[size(f,1)]
    "Absolute value/magnitude vector (either without unit or in dB)";
  input SI.Conversions.NonSIunits.Angle_deg phi[size(f,1)] "Angles in degree";

  input Boolean autoRange=true
    "= true, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min=0.1
    "Minimum frequency value, if autoRange = false";
  input Modelica.SIunits.Frequency f_max=10
    "Maximum frequency value, if autoRange = false";

  input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
  input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

  input Boolean Hz=true
    "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
  input Boolean dB=false
    "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));

  input Plot.Records.Diagram diagram "Diagram layout" annotation(Dialog);
  input Plot.Records.Device device=Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);
protected
  Boolean OK;
  Integer window=0;
  SI.Angle phi_old;
  Plot.Records.Curve curves[2];
  Integer i;
  Plot.Records.Diagram diagram2[2];
algorithm
  // Plot computed frequency response
  diagram2 := fill(diagram, 2);
  i := 0;
  if magnitude then
    i := i + 1;
    curves[i] := Plot.Records.Curve(
      x=f,
      y=a,
      autoLine=true);
    diagram2[i].curve := {curves[i]};
    diagram2[i].yLabel := if dB then "magnitude [dB]" else "magnitude";
    if phase then
       diagram2[i].xLabel:="";
    end if;
    if dB then
       diagram2[i].logY := false;
    end if;
  end if;

  if phase then
    i := i + 1;
    curves[i] := Plot.Records.Curve(
      x=f,
      y=phi,
      autoLine=true);
    diagram2[i].curve := {curves[i]};
    diagram2[i].yLabel := "phase [deg]";
    diagram2[i].logY := false;
    if magnitude then
      diagram2[i].heading:="";
    end if;
  end if;

  if not Hz then
     diagram2[i].xLabel:="Angular frequency [rad/s]";
  end if;

  if magnitude and phase then
    Plot.diagramVector(diagram2, device);
  else
    Plot.diagram(diagram2[1], device);
  end if;

  annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>bode</b>(zp)
   or
ZerosAndPoles.Plot.<b>bode</b>(
  zp,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  diagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the bode-diagram of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(p^2 + 5*p + 7)/(p + 2)/(p + 3);

<b>algorithm</b>
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(zp)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</blockquote>
</html>"));
end frequencyResponsePlot;
