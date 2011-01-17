within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Example_UKF
  import Modelica_LinearSystems2;

  Modelica.Blocks.Sources.Sine sine(
    startTime=1,
    amplitude=0.1,
    freqHz=2,
    offset=0) "Force on grap"
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica_LinearSystems2.Controller.Noise noise(
    y_min=-0.1,
    y_max=0.1,
    sampleFactor=3)
    annotation (Placement(transformation(extent={{-90,-90},{-70,-70}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete,
    initType=Modelica_LinearSystems2.Controller.Types.Init.InitialState,
    sampleTime=0.001)
             annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica_LinearSystems2.Controller.Noise noise1(
    y_min=-0.003,
    y_max=0.003,
    sampleFactor=1)
    annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
  Modelica_LinearSystems2.WorkInProgress.MPC.Examples.SystemModel systemModel(
      system(
      A=[-1.2822,0,0.98,0; 0,0,1,0; -5.4293,0,-1.8366,0; -128.2,128.2,0,0],
      B=[-0.3; 0; -17; 0],
      C=[0,1,0,0; 0,0,0,1; -128.2,128.2,0,0],
      D=[0; 0; 0]))
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Routing.Multiplex3 multiplex3_1
    annotation (Placement(transformation(extent={{20,-40},{40,-20}})));

  Modelica_LinearSystems2.Controller.UKF UKF(
    beta=2,
    Q=1e-5*identity(size(UKF.x_est_init, 1)),
    G=identity(size(UKF.x_est_init, 1)),
    kappa=1,
    alpha=0.1,
    P_init=0.5*identity(size(UKF.x_est_init, 1)),
    redeclare function F_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.fSigmaLinear,
    redeclare function H_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.hSigmaLinear,
    x_est_init={0,0,0,0},
    R=diagonal({0.1,1,2.13e-1}),
    sampleFactor=10)
    annotation (Placement(transformation(extent={{60,-20},{80,0}})));
equation
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-50,-76},{-50,-80},{-69,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-22,-36},{-50,-36},{-50,-40},{-69,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(systemModel.u[1], sine.y) annotation (Line(
      points={{-62,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(systemModel.y[1], add.u1) annotation (Line(
      points={{-39,10},{-28,10},{-28,-24},{-22,-24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u1, systemModel.y[2]) annotation (Line(
      points={{-22,-64},{-28,-64},{-28,10},{-39,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, multiplex3_1.u1[1]) annotation (Line(
      points={{1,-30},{6,-30},{6,-23},{18,-23}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.u2[1], add1.y) annotation (Line(
      points={{18,-30},{10,-30},{10,-70},{1,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.u3[1], systemModel.y[3]) annotation (Line(
      points={{18,-37},{10,-37},{10,10},{-39,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.y, UKF.y_measure) annotation (Line(
      points={{41,-30},{50,-30},{50,-15},{58,-15}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, UKF.u[1]) annotation (Line(
      points={{-79,10},{-70,10},{-70,34},{30,34},{30,-5},{58,-5}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=5),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResultsEKF.mos"
        "plotSetup"));
end Example_UKF;
