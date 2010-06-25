within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_MBD_UKF
  import Modelica_LinearSystems2;

  Modelica_LinearSystems2.Controller.UKF UKF(
    beta=2,
    x_est_init={0.1,0.1,1,0.2},
    R=diagonal({0.1,2.13e-1}),
    Q=1e-5*identity(size(UKF.x_est_init, 1)),
    G=identity(size(UKF.x_est_init, 1)),
    kappa=1,
    redeclare function F_function = 
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.fSigma,
    alpha=0.1,
    redeclare function H_function = 
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.hSigma,
    P_init=0.5*identity(size(UKF.x_est_init, 1))) 
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  CraneMultibody craneWithEquations3_1(d=100, J=0) 
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Sources.Sine sine(
    amplitude=2000,
    freqHz=0.5,
    offset=-20,
    startTime=1) "Force on grap" 
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica_LinearSystems2.Controller.Noise noise(
    sampleFactor=3,
    y_min=-0.1,
    y_max=0.1) 
    annotation (Placement(transformation(extent={{-90,-80},{-70,-60}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1 
    annotation (Placement(transformation(extent={{20,-30},{40,-50}})));
  Modelica.Blocks.Math.Add add 
    annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
  Modelica.Blocks.Math.Add add1 
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete, sampleTime=
       0.005) 
             annotation (Placement(transformation(extent={{80,80},{100,100}})));
  Modelica_LinearSystems2.Controller.Noise noise1(
    sampleFactor=8,
    y_min=-0.03,
    y_max=0.03) 
    annotation (Placement(transformation(extent={{-90,-40},{-70,-20}})));
equation
  connect(sine.y, craneWithEquations3_1.force) annotation (Line(
      points={{-79,10},{-62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(UKF.u[1], sine.y)            annotation (Line(
      points={{58,5},{46,5},{46,40},{-70,40},{-70,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-44,-76},{-44,-70},{-69,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(craneWithEquations3_1.y1, add1.u1) annotation (Line(
      points={{-39,1},{-30,1},{-30,-64},{-22,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-22,-16},{-46,-16},{-46,-30},{-69,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, craneWithEquations3_1.x1) annotation (Line(
      points={{-22,-4},{-26,-4},{-26,19},{-39,19}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2[1], add.y) annotation (Line(
      points={{18,-34},{12,-34},{12,-10},{1,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u1[1], add1.y) annotation (Line(
      points={{18,-46},{10,-46},{10,-70},{1,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.y, UKF.y_measure) annotation (Line(
      points={{41,-40},{46,-40},{46,-5},{58,-5}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=20),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResultsUKF.mos"
        "plotResults"));
end Pendular_Testbench_MBD_UKF;
