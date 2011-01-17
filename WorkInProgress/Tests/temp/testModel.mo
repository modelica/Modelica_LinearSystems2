within Modelica_LinearSystems2.WorkInProgress.Tests.temp;
model testModel
  blockFunction blockFunction1
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Modelica.Blocks.Sources.Step step
    annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
  Modelica.Blocks.Sources.Sine sine
    annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1
    annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_2
    annotation (Placement(transformation(extent={{-60,10},{-40,-10}})));
equation
  connect(step.y, multiplex2_1.u1[1]) annotation (Line(
      points={{-79,30},{-70,30},{-70,26},{-62,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2[1], sine.y) annotation (Line(
      points={{-62,14},{-70,14},{-70,-10},{-79,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_2.u2[1], step.y) annotation (Line(
      points={{-62,6},{-72,6},{-72,30},{-79,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_2.u1[1], sine.y) annotation (Line(
      points={{-62,-6},{-72,-6},{-72,-10},{-79,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.y, blockFunction1.u1) annotation (Line(
      points={{-39,20},{-22,20},{-22,14},{-2,14}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(blockFunction1.u2, multiplex2_2.y) annotation (Line(
      points={{-2,4},{-22,4},{-22,0},{-39,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end testModel;
