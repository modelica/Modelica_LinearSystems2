within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_MBD_USKF

  SR_UnscentedKalman2 UKF(
    beta=2,
    Ts=0.005,
    alpha=0.01,
    xm_start={0.1,0.1,1,0.2}) 
    annotation (Placement(transformation(extent={{60,2},{80,22}})));
  CraneMultibody craneWithEquations3_1(d=100, J=0) 
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Sources.Sine sine(amplitude=1000, freqHz=0.1,
    offset=-50) "Force on grap" 
    annotation (Placement(transformation(extent={{-98,0},{-78,20}})));
  Modelica_LinearSystems2.Controller.Noise noise(y_min=-0.05, y_max=0.05) 
    annotation (Placement(transformation(extent={{-92,-70},{-72,-50}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1 
    annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
  Modelica.Blocks.Math.Add add 
    annotation (Placement(transformation(extent={{-18,-40},{2,-20}})));
  Modelica.Blocks.Math.Add add1 
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete, sampleTime=
       0.005) 
             annotation (Placement(transformation(extent={{80,78},{100,98}})));
  Modelica_LinearSystems2.Controller.Noise noise1(y_min=-0.0001, y_max=0.0001) 
    annotation (Placement(transformation(extent={{-90,-34},{-70,-14}})));
equation
  connect(sine.y, craneWithEquations3_1.force) annotation (Line(
      points={{-77,10},{-62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(UKF.u[1], sine.y)            annotation (Line(
      points={{58,16},{16,16},{16,42},{-70,42},{-70,10},{-77,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-58,-76},{-58,-60},{-71,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(UKF.y[:, 1], multiplex2_1.y)            annotation (Line(
      points={{58,8},{48,8},{48,-40},{41,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(craneWithEquations3_1.y1, add1.u1) annotation (Line(
      points={{-39,1},{-30,1},{-30,-64},{-22,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-20,-36},{-50,-36},{-50,-24},{-69,-24}},
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
  annotation (Diagram(graphics),
    experiment(StopTime=50),
    experimentSetupOutput,
    Commands(file="WorkInProgress/Tests/Examples/plotResultsUKF.mos"
        "plotResultsUKF"));
end Pendular_Testbench_MBD_USKF;
