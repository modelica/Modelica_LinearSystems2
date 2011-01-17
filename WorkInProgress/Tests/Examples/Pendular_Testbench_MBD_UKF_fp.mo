within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_MBD_UKF_fp
  import Modelica_LinearSystems2;

  CraneMultibody craneWithEquations3_1(     J=0, d=100)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Sources.Sine sine(
    amplitude=2000,
    freqHz=0.5,
    offset=-20,
    startTime=1) "Force on grap"
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica_LinearSystems2.Controller.Noise noise(
    y_min=-0.1,
    y_max=0.1,
    sampleFactor=3)
    annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1
    annotation (Placement(transformation(extent={{20,-40},{40,-60}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete, sampleTime=
       0.005)
             annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica_LinearSystems2.Controller.Noise noise1(
    y_min=-0.03,
    y_max=0.03,
    sampleFactor=8)
    annotation (Placement(transformation(extent={{-90,-46},{-70,-26}})));
  Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter.EKF EKF(
    Q=1e-5*identity(size(EKF.x_est_init, 1)),
    G=identity(size(EKF.x_est_init, 1)),
    x_est_init={0.1,0.1,1,0.2},
    R=diagonal({0.1,2.13*1e-1}),
    M_init=0.5*identity(size(EKF.x_est_init, 1)),
    sampleFactor=1,
    redeclare function ekfFunction =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.ekfSystem_pendular)
    annotation (Placement(transformation(extent={{60,-60},{80,-40}})));

  Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter.UKF uKF(
    redeclare function F_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.fSigma,
    redeclare function H_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.hSigma,
    x_est_init={0.1,0.1,1,0.2},
    Q=1e-5*identity(size(EKF.x_est_init, 1)),
    R=diagonal({0.1,2.13*1e-1}),
    G=identity(size(EKF.x_est_init, 1)),
    P_init=0.5*identity(size(EKF.x_est_init, 1)))
    annotation (Placement(transformation(extent={{60,-20},{80,0}})));
  Modelica_LinearSystems2.WorkInProgress.Controller.KalmanFilter.UKF_SR uKF_SR(
    redeclare function F_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.fSigma,
    redeclare function H_function =
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.hSigma,
    x_est_init={0.1,0.1,1,0.2},
    Q=1e-5*identity(size(EKF.x_est_init, 1)),
    R=diagonal({0.1,2.13*1e-1}),
    G=identity(size(EKF.x_est_init, 1)),
    P_init=0.5*identity(size(EKF.x_est_init, 1)))
    annotation (Placement(transformation(extent={{60,20},{80,40}})));
equation
  connect(sine.y, craneWithEquations3_1.force) annotation (Line(
      points={{-79,10},{-62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-69,-76}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(craneWithEquations3_1.y1, add1.u1) annotation (Line(
      points={{-39,1},{-30,1},{-30,-64},{-22,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-22,-36},{-69,-36}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, craneWithEquations3_1.x1) annotation (Line(
      points={{-22,-24},{-26,-24},{-26,19},{-39,19}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2[1], add.y) annotation (Line(
      points={{18,-44},{12,-44},{12,-30},{1,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u1[1], add1.y) annotation (Line(
      points={{18,-56},{10,-56},{10,-70},{1,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(EKF.u[1], sine.y) annotation (Line(
      points={{58,-45},{50,-45},{50,40},{-68,40},{-68,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(EKF.y_measure, multiplex2_1.y) annotation (Line(
      points={{58,-55},{50,-55},{50,-50},{41,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(uKF.u[1], sine.y) annotation (Line(
      points={{58,-5},{40,-5},{40,16},{18,16},{18,40},{-68,40},{-68,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(uKF.y_measure, multiplex2_1.y) annotation (Line(
      points={{58,-15},{42,-15},{42,-50},{41,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(uKF_SR.u[1], sine.y) annotation (Line(
      points={{58,35},{-70,35},{-70,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(uKF_SR.y_measure, multiplex2_1.y) annotation (Line(
      points={{58,25},{41,25},{41,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=50),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResultsEKF.mos"
        "plotSetup"));
end Pendular_Testbench_MBD_UKF_fp;
