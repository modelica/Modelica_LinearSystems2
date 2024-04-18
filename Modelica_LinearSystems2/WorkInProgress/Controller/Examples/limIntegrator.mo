within Modelica_LinearSystems2.WorkInProgress.Controller.Examples;
model limIntegrator "linIntegrator"
  extends Modelica.Icons.Example;

  import Modelica_LinearSystems2;
  Modelica_LinearSystems2.WorkInProgress.Controller.LimIntegrator limIntegrator(
    y_start=1.5,
    limitsAtInit=false,
    withDelay=true,
    initType=Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault.InitialState)
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Modelica.Blocks.Sources.Sine sine(f=1, amplitude=5)
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
  Modelica.Blocks.Sources.Pulse pulse(
    period=0.25,
    offset=1.3,
    amplitude=1)
    annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
  Modelica.Blocks.Sources.Pulse pulse1(
    amplitude=0.1,
    offset=0.5,
    period=0.34)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  inner Modelica_LinearSystems2.Controllers.SampleClock sampleClock(
    initType=Modelica_LinearSystems2.Controllers.Types.Init.InitialState,
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockType.Continuous,
    sampleTime=0.02)
             annotation (Placement(transformation(extent={{54,50},{74,70}})));

  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-8,-74},{12,-54}})));
equation
  connect(sine.y, limIntegrator.u) annotation (Line(
      points={{-59,10},{-2,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limIntegrator.limit1, pulse.y) annotation (Line(
      points={{-2,18},{-20,18},{-20,50},{-59,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limIntegrator.limit2, pulse1.y) annotation (Line(
      points={{-2,2},{-14,2},{-14,-30},{-59,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integratorXX.u, sine.y) annotation (Line(
      points={{22,-38},{-22,-38},{-22,10},{-59,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integratorXX.lowerLimit, pulse1.y) annotation (Line(
      points={{22,-46},{-22,-46},{-22,-30},{-59,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integratorXX.upperLimit, pulse.y) annotation (Line(
      points={{22,-30},{-20,-30},{-20,50},{-59,50}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    __Dymola_experimentSetupOutput);
end limIntegrator;
