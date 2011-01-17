within Modelica_LinearSystems2.WorkInProgress.Controller.Examples;
model limIntegrator "linIntegrator"
  import Modelica_LinearSystems2;
  Modelica_LinearSystems2.WorkInProgress.Controller.LimIntegrator limIntegrator(
    y_start=1.5,
    limitsAtInit=false,
    withDelay=true,
    initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.InitialState)
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Modelica.Blocks.Sources.Sine sine(freqHz=1, amplitude=5)
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
  Modelica.Blocks.Sources.Pulse pulse(
    period=1,
    amplitude=0.1,
    offset=1.5)
    annotation (Placement(transformation(extent={{-76,40},{-56,60}})));
  Modelica.Blocks.Sources.Pulse pulse1(
    period=1.3,
    amplitude=0.1,
    offset=0.5)
    annotation (Placement(transformation(extent={{-76,-48},{-56,-28}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(
    blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Discrete,
    sampleTime=0.02,
    initType=Modelica_LinearSystems2.Controller.Types.Init.InitialState)
             annotation (Placement(transformation(extent={{54,50},{74,70}})));
  Modelica_LinearSystems2.Controller.LimPID PID(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    kp=0,
    yMin=0.5,
    y_start=1.5,
    initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.NoInit,

    yMax=1) annotation (Placement(transformation(extent={{28,-42},{48,-22}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-8,-74},{12,-54}})));
equation
  connect(sine.y, limIntegrator.u) annotation (Line(
      points={{-59,10},{-2,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limIntegrator.limit1, pulse.y) annotation (Line(
      points={{-2,18},{-24,18},{-24,50},{-55,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limIntegrator.limit2, pulse1.y) annotation (Line(
      points={{-2,2},{-14,2},{-14,-38},{-55,-38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID.u_s, sine.y) annotation (Line(
      points={{26,-32},{-18,-32},{-18,10},{-59,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, PID.u_m) annotation (Line(
      points={{13,-64},{38,-64},{38,-44}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(graphics),
    experiment(StopTime=5),
    __Dymola_experimentSetupOutput);
end limIntegrator;
