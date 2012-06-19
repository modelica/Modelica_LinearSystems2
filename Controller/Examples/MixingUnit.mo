within Modelica_LinearSystems2.Controller.Examples;
model MixingUnit
  "Example of system control with inverse model in a controller with two degrees of freedom"
  import Modelica_LinearSystems2;
  extends Modelica.Icons.Example;
  extends Templates.TwoDOFinverseModelController(
    redeclare Examples.Components.MixingUnit plant_inv(
        mixingUnit(
        c(start=c_start, fixed=true),
        T_c(start=T_c_start, fixed=true),
        T(start=T_start, fixed=false))),
    redeclare Examples.Components.MixingUnit plant(
        mixingUnit(c(start=c_start, fixed=true), T(start=T_start, fixed=true))),
    filter(
      order=3,
      normalized=false,
      f_cut=freq,
      initType=Types.InitWithGlobalDefault.NoInit),
    redeclare Controller.PI controller(k=10, T=10));

  import SI = Modelica.SIunits;
  parameter Real x10 = 0.42;
  parameter Real x10_inv = 0.6;
  parameter Real x20 = 0.01;
  parameter Real u0 = -0.0224;
  parameter SI.Frequency freq = 1/300 "Critical frequency of filter";
  final parameter Real c0 = 0.848;
  final parameter SI.Temperature T0 = 308.5;
  final parameter Real c_start(unit="mol/l") = c0*(1-x10);
  final parameter Real c_inv_start(unit="mol/l") = c0*(1-x10_inv);
  final parameter SI.Temperature T_start = T0*(1+x20);
  final parameter Real c_high_start(unit="mol/l") = c0*(1-0.72);
  final parameter Real T_c_start = T0*(1+u0);

  Modelica.Blocks.Sources.Step step1(height=c_high_start - c_start, offset=
        c_start,
    startTime=25)
    annotation (Placement(transformation(extent={{-134,10},{-114,30}},
          rotation=0)));
  inner Controller.SampleClock sampleClock
    annotation (Placement(transformation(extent={{80,60},{100,80}})));
equation
  connect(step1.y, filter.u) annotation (Line(
      points={{-113,20},{-102,20}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=500),
    Commands(
      file="modelica://Modelica_LinearSystems2/Resources/Scripts/Dymola/Controllers/Examples/MixingUnit_plot.mos"
        "Plot Results",
      file(
        ensureSimulated=true,
        partOfCheck=true)=
        "modelica://Modelica_LinearSystems2/Resources/Scripts/Dymola/Controllers/Examples/MixingUnit_plot.mos"
        "Simulate and Plot Results"),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-140,-100},{140,
            100}}), graphics),
    Documentation(info="<html>
<p>This example demonstrates the usage of the control structure template <i>Modelica_Controller.Templates.TwoDOFinverseModelController2</i> to control a system by using of a inverse system model in the forward path. The controlled system is a mixing unit described in [1]. </p>
<p>Within Dymola simulation tool the &QUOT;Commands / Simulate and Plot Results&QUOT; selection yields the following simulation result. </p>
<h4>Simulation results </h4>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/Examples/MixingUnit_results.png\"/> </p>
<h4><a name=\"References\">R</a>eferences</h4>
<pre>  [1] F&ouml;llinger O., &QUOT;Nichtlineare Regelungen I&QUOT;,
      Oldenbourg Verlag, 8. Auflage.</pre>
</html>"));
end MixingUnit;
