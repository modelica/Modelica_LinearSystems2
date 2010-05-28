within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Pendular_Testbench_WEQ

  ExtendedKalman extendedKalman(Ts=0.001) 
    annotation (Placement(transformation(extent={{60,0},{80,20}})));
  CraneWithEquations craneWithEquations3_1(phi(
      displayUnit="rad",
      fixed=true,
      start=0), m_load=3000) 
    annotation (Placement(transformation(extent={{-56,0},{-36,20}})));
  Modelica.Blocks.Sources.Sine sine(amplitude=5000, freqHz=0.5) "Force on grap"
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
      points={{-77,10},{-58,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(extendedKalman.u[1], sine.y) annotation (Line(
      points={{58,14},{16,14},{16,42},{-70,42},{-70,10},{-77,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u2, noise.y) annotation (Line(
      points={{-22,-76},{-58,-76},{-58,-60},{-71,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(extendedKalman.y[:, 1], multiplex2_1.y) annotation (Line(
      points={{58,6},{48,6},{48,-40},{41,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(craneWithEquations3_1.y1, add1.u1) annotation (Line(
      points={{-35,16},{-30,16},{-30,-64},{-22,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, noise1.y) annotation (Line(
      points={{-20,-36},{-50,-36},{-50,-24},{-69,-24}},
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
  connect(craneWithEquations3_1.y2, add.u1) annotation (Line(
      points={{-35,4},{-34,4},{-34,-24},{-20,-24}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=10),
    experimentSetupOutput);
end Pendular_Testbench_WEQ;
