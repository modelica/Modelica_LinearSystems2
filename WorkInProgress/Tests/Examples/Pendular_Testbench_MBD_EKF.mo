within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_MBD_EKF
  import Modelica_LinearSystems2;

  Modelica_LinearSystems2.Controller.EKF EKF(
    redeclare function ekfFunction = 
        Modelica_LinearSystems2.WorkInProgress.Tests.Examples.ekfSystem_pendular,
    nu=1,
    Q=1e-5*identity(size(EKF.x_est_init, 1)),
    G=identity(size(EKF.x_est_init, 1)),
    x_est_init={0.1,0.1,1,0.2},
    R=diagonal({0.1,2.13*1e-1}),
    M_init=0.5*identity(size(EKF.x_est_init, 1))) 
    annotation (Placement(transformation(extent={{60,-20},{80,0}})));

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
    annotation (Placement(transformation(extent={{-90,-90},{-70,-70}})));
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
    annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
equation
  connect(sine.y, craneWithEquations3_1.force) annotation (Line(
      points={{-79,10},{-62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(EKF.u[1], sine.y)            annotation (Line(
      points={{58,-5},{20,-5},{20,40},{-70,40},{-70,10},{-79,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-50,-76},{-50,-80},{-69,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(craneWithEquations3_1.y1, add1.u1) annotation (Line(
      points={{-39,1},{-30,1},{-30,-64},{-22,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-22,-36},{-50,-36},{-50,-40},{-69,-40}},
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
  connect(EKF.y_measure, multiplex2_1.y)            annotation (Line(
      points={{58,-15},{49,-15},{49,-50},{41,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=50),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResultsEKF.mos"
        "plotSetup"));
end Pendular_Testbench_MBD_EKF;
