within Modelica_LinearSystems2.Controller.Examples;
model MixingUnit
  "Example of system control with inverse model in a controller with two degrees of freedom"
  extends Modelica.Icons.Example;
  extends Templates.TwoDOFinverseModelController(
    redeclare Examples.Components.MixingUnit plant_inv(mixingUnit(
        c(start=c_start, fixed=true),
        T_c(start=T_c_start, fixed=true),
        T(start=T_start, fixed=true))),
    redeclare Examples.Components.MixingUnit plant(
        mixingUnit(c(start=c_start, fixed=true), T(start=T_start, fixed=true))),
    filter(
      order=3,
      normalized=false,
      f_cut=freq,
      initType=Types.InitWithGlobalDefault.NoInit),
    redeclare Controller.PI controller(k=10, T=10,
      initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.InitialState));

  import SI = Modelica.SIunits;
  parameter Real x10 = 0.42
    "Initial value of state x1 (related concentration of substance A in tank)";
  parameter Real x10_inv = 0.6 "Initial value of state x1 of inverted model";
  parameter Real x20 = 0.01
    "Initial value of state x2 (related temperature in tank)";
  parameter Real u0 = -0.0224
    "Initial related temperature of cooling medium [-]";
  parameter SI.Frequency freq = 1/300 "Critical frequency of filter";

  final parameter Real c0 = 0.848
    "Nominal concentration of substance A on intake";
  final parameter SI.Temperature T0 = 308.5
    "Nominal temperature of substance A on intake";
  final parameter Real c_start(unit="mol/l") = c0*(1-x10)
    "Initial concentration of substance A in tank";
  final parameter Real c_inv_start(unit="mol/l") = c0*(1-x10_inv)
    "Initial concentration of substance A in tank";
  final parameter SI.Temperature T_start = T0*(1+x20)
    "Initial temperature in tank";
  final parameter Real c_high_start(unit="mol/l") = c0*(1-0.72)
    "Concentration change height";
  final parameter SI.Temperature T_c_start = T0*(1+u0)
    "Initial temperature of cooling medium";

  Modelica.Blocks.Sources.Step step1(
    height=c_high_start - c_start,
    offset=c_start,
    startTime=25)
    annotation (Placement(transformation(extent={{-120,10},{-100,30}},
          rotation=0)));
  inner Controller.SampleClock sampleClock
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  connect(step1.y, filter.u) annotation (Line(
      points={{-99,20},{-92,20}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=500),
    __Dymola_Commands(
      file="modelica://Modelica_LinearSystems2/Resources/Scripts/Dymola/Controllers/Examples/MixingUnit_plot.mos"
        "Plot Results",
      file(
        ensureSimulated=true,
        partOfCheck=true)=
        "modelica://Modelica_LinearSystems2/Resources/Scripts/Dymola/Controllers/Examples/MixingUnit_plot.mos"
        "Simulate and Plot Results"),
    Documentation(info="<html>
<p>
This example demonstrates the usage of the control structure template
<i>Modelica_Controller.Templates.TwoDOFinverseModelController2</i> to
control a system by using of an inverse system model in the forward path.
The controlled system is a mixing unit described in [1].
See also model of
<a href=\"modelica://Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit1\">MixingUnit1</a>
for more details.
</p>
<p>
Within Dymola simulation tool the &quot;Commands / Simulate and Plot Results&quot;
selection yields the following simulation result.
</p>

<h4>Simulation results </h4>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/Examples/MixingUnit_results.png\"/> </p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] F&ouml;llinger O. (1998):</dt>
<dd> <b>Nichtlineare Regelungen I</b>.
     8th Edition, Oldenbourg Verlag M&uuml;nchen.<br>&nbsp;</dd>
</dl>
</html>"),
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-100},{120,
            100}}), graphics));
end MixingUnit;
