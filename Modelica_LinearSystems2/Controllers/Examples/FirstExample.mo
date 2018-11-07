within Modelica_LinearSystems2.Controllers.Examples;
model FirstExample "First example to demonstrate representative block"
  extends Modelica.Icons.Example;
  import Modelica_LinearSystems2;

  parameter Modelica.SIunits.AngularFrequency w=10 "Undamped natural frequency";
  parameter Real D=0.1 "Damping ratio";

  Modelica.Blocks.Sources.Step step(
    startTime=0.5,
    height=1.2,
    offset=0.2)                      annotation (
      Placement(transformation(extent={{-80,0},{-60,20}})));
  Modelica_LinearSystems2.Controllers.StateSpace stateSpace(
    x_start={0.1,0},
    initType=Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault.InitialState,
    system=Modelica_LinearSystems2.StateSpace(
      A=[0,1; -w*w,-2*w*D],
      B=[0; w*w],
      C=[1,0],
      D=[0]),
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
      annotation(Placement(transformation(extent={{-20,40},{0,60}})));

  TransferFunction transferFunction(
    system(n={1,2}, d={1,2,3}),
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
    annotation (Placement(transformation(extent={{-20,0},{0,20}})));
  ZerosAndPoles zerosAndPoles(
    system(
      n1={1},
      n2=fill(0, 0, 2),
      d1=fill(0, 0),
      d2=[1,1; 1,1]),
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
  inner SampleClock sampleClock(
    sampleTime=0.1,
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockType.Continuous)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
equation
  connect(step.y, stateSpace.u[1])      annotation (Line(
      points={{-59,10},{-40,10},{-40,50},{-22,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(transferFunction.u, step.y) annotation (Line(
      points={{-22,10},{-59,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.u, step.y) annotation (Line(
      points={{-22,-30},{-40,-30},{-40,10},{-59,10}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (    experiment(StopTime=5));
end FirstExample;
