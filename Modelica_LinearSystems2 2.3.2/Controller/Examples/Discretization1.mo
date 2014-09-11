within Modelica_LinearSystems2.Controller.Examples;
model Discretization1 "Demonstrates the discretization methods"
  extends Modelica.Icons.Example;

  parameter Real w=20 "Angular frequency";
  parameter Real D=0.1 "Damping";

  Modelica_LinearSystems2.Controller.SecondOrder continuous(w=w, D=D)
    annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(sampleTime=
        0.01)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica.Blocks.Sources.Step step(
    height=1.2,
    offset=0.2,
    startTime=0.1)                   annotation (
      Placement(transformation(extent={{-80,40},{-60,60}})));
  Modelica_LinearSystems2.Controller.SecondOrder explicitEuler(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ExplicitEuler)
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

  Modelica_LinearSystems2.Controller.SecondOrder implicitEuler(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ImplicitEuler)
    annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));

  Modelica_LinearSystems2.Controller.SecondOrder trapezoid(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.Trapezoidal)
    annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));

  Modelica_LinearSystems2.Controller.SecondOrder impulseExact(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.ImpulseExact)
    annotation (Placement(transformation(extent={{20,20},{40,40}})));

  Modelica_LinearSystems2.Controller.SecondOrder stepExact(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.StepExact)
    annotation (Placement(transformation(extent={{20,-20},{40,0}})));

  Modelica_LinearSystems2.Controller.SecondOrder rampExact(
    w=w,
    D=D,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    methodType=Modelica_LinearSystems2.Controller.Types.MethodWithGlobalDefault.RampExact)
    annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

equation
  connect(step.y, continuous.u) annotation (Line(
      points={{-59,50},{-42,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, explicitEuler.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,10},{-42,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, implicitEuler.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,-30},{-42,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, trapezoid.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,-70},{-42,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, impulseExact.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,80},{0,80},{0,30},{18,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, stepExact.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,80},{0,80},{0,-10},{18,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, rampExact.u) annotation (Line(
      points={{-59,50},{-52,50},{-52,80},{0,80},{0,-50},{18,-50}},
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
end Discretization1;
