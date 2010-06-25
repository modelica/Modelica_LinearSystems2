within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_MBD_SR_UKF
  import Modelica_LinearSystems2;

  Modelica_LinearSystems2.Controller.UKF_SR UKF_SR(
   beta=2,
    x_est_init={0.1,0.1,1,0.2},
    R=diagonal({0.1,2.13e-1}),
    Q=1e-5*identity(size(UKF_SR.x_est_init, 1)),
    G=identity(size(UKF_SR.x_est_init, 1)),
    redeclare function F_function = 
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.fSigma,
    alpha=0.1,
    redeclare function H_function = 
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.hSigma,
    P_init=0.5*identity(size(UKF_SR.x_est_init, 1)),
    kappa=0) 
    annotation (Placement(transformation(extent={{60,0},{80,20}})));
  CraneMultibody craneWithEquations3_1(d=100, J=0) 
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Sources.Sine sine(
    freqHz=0.5,
    amplitude=2000,
    offset=-20,
    startTime=1) "Force on grap" 
    annotation (Placement(transformation(extent={{-98,0},{-78,20}})));
  Modelica_LinearSystems2.Controller.Noise noise(y_min=-0.1, y_max=0.1,
    sampleFactor=3) 
    annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1 
    annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
  Modelica.Blocks.Math.Add add 
    annotation (Placement(transformation(extent={{-18,-40},{2,-20}})));
  Modelica.Blocks.Math.Add add1 
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete, sampleTime=
       0.005) 
             annotation (Placement(transformation(extent={{80,80},{100,100}})));
  Modelica_LinearSystems2.Controller.Noise noise1(
    y_min=-0.03,
    y_max=0.03,
    sampleFactor=8) 
    annotation (Placement(transformation(extent={{-90,-46},{-70,-26}})));
equation
  connect(sine.y, craneWithEquations3_1.force) annotation (Line(
      points={{-77,10},{-62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(UKF_SR.u[1], sine.y)         annotation (Line(
      points={{58,15},{20,15},{20,40},{-70,40},{-70,10},{-77,10}},
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
      points={{-20,-36},{-69,-36}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, craneWithEquations3_1.x1) annotation (Line(
      points={{-20,-24},{-26,-24},{-26,19},{-39,19}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2[1], add.y) annotation (Line(
      points={{18,-46},{12,-46},{12,-30},{3,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u1[1], add1.y) annotation (Line(
      points={{18,-34},{10,-34},{10,-70},{1,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(UKF_SR.y_measure, multiplex2_1.y) annotation (Line(
      points={{58,5},{48,5},{48,-40},{41,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=50),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResults_SR_UKF.mos"
        "plotResults"));
end Pendular_Testbench_MBD_SR_UKF;
