within Modelica_LinearSystems2.Controllers.Examples;
model Discretization2
  "Demonstrates the discretization method for impuls exact discretization"
  extends Modelica.Icons.Example;

  parameter Real w=20 "Angular frequency";
  parameter Real D=0.1 "Damping";

  inner Modelica_LinearSystems2.Controllers.SampleClock sampleClock(
    sampleTime=0.01)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));

  Modelica_LinearSystems2.Controllers.SecondOrder impulseExact(
    D=D,
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controllers.Types.MethodWithGlobalDefault.ImpulseExact,
    w=w)
    annotation (Placement(transformation(extent={{0,-30},{20,-10}})));

  Modelica_LinearSystems2.Controllers.SecondOrder continuous(D=D, w=w)
    annotation (Placement(transformation(extent={{0,10},{20,30}})));
  Derivative derivative(T=1e-8)
    annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
  Modelica.Blocks.Sources.Pulse pulse(
    startTime=0.1,
    period=1,
    width=sampleClock.sampleTime*100)
    annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
  Modelica.Blocks.Sources.Step step1(
    startTime=0.1,
    height=1,
    offset=0) annotation (Placement(transformation(extent={{-80,10},{-60,30}})));

equation
  connect(pulse.y, impulseExact.u) annotation (Line(
      points={{-59,-20},{-2,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(continuous.u, derivative.y) annotation (Line(
      points={{-2,20},{-19,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step1.y, derivative.u) annotation (Line(
      points={{-59,20},{-42,20}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(Tolerance=1e-006),
    Documentation(info="<html>
<p>
Demonstrates the different discretization methods by simulating the step
response of a second order system as continuous system and as discrete system
with the supported discretization methods. The step starts with an offset at 0.1 s
to demonstrate the steady-state initialization.
</p>
</html>"));
end Discretization2;
