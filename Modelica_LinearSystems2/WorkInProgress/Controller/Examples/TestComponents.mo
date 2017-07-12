within Modelica_LinearSystems2.WorkInProgress.Controller.Examples;
model TestComponents "test all Controller blocks"
  extends Modelica.Icons.Example;

  parameter Real w=10;
  parameter Real D=0.1;
  Modelica.Blocks.Sources.Step step(
    startTime=0.5,
    height=1.2,
    offset=0.2)                      annotation (extent=[-80,40; -60,60],
      Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica_LinearSystems2.Controllers.StateSpace stateSpace(
    x_start={0.1,0},
    initType=Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault.InitialState,
    system(
      A=[0,1; -w*w,-2*w*D],
      B=[0; w*w],
      C=[1,0],
      D=[0]),
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
                                                              annotation (extent=[-40,40;
        -20,60], Placement(transformation(extent={{-20,220},{0,240}})));

  Modelica_LinearSystems2.Controllers.TransferFunction transferFunction(system(n=
         {1,2}, d={1,2,3}), blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
    annotation (Placement(transformation(extent={{-20,190},{0,210}})));
  Modelica_LinearSystems2.Controllers.ZerosAndPoles zerosAndPoles(system(
      n1={1},
      n2=fill(
          0,
          0,
          2),
      d1=fill(0, 0),
      d2=[1,1; 1,1]), blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption)
    annotation (Placement(transformation(extent={{-20,160},{0,180}})));
  inner Modelica_LinearSystems2.Controllers.SampleClock sampleClock(
    sampleTime=0.1,
    initType=Modelica_LinearSystems2.Controllers.Types.Init.InitialState,
    blockType=Modelica_LinearSystems2.Controllers.Types.BlockType.Discrete)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica_LinearSystems2.Controllers.Filter filter
    annotation (Placement(transformation(extent={{-20,130},{0,150}})));
  Modelica_LinearSystems2.Controllers.FilterFIR filter1
    annotation (Placement(transformation(extent={{-20,100},{0,120}})));
  Modelica_LinearSystems2.Controllers.Integrator integrator
    annotation (Placement(transformation(extent={{-20,70},{0,90}})));
  Modelica_LinearSystems2.Controllers.Derivative derivative
    annotation (Placement(transformation(extent={{-20,40},{0,60}})));
  Modelica_LinearSystems2.Controllers.FirstOrder firstOrder
    annotation (Placement(transformation(extent={{-20,10},{0,30}})));
  Modelica_LinearSystems2.Controllers.SecondOrder secondOrder
    annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
  Modelica_LinearSystems2.Controllers.PI pI
    annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
  Modelica_LinearSystems2.Controllers.PID pID(initType=Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault.InitialOutput)
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  Modelica_LinearSystems2.Controllers.LimPID PID(initType=
        Modelica_LinearSystems2.Controllers.Types.InitWithGlobalDefault.NoInit)
    annotation (Placement(transformation(extent={{-20,-110},{0,-90}})));
  Modelica_LinearSystems2.Controllers.UnitDelay unitDelay
    annotation (Placement(transformation(extent={{-20,-160},{0,-140}})));
  Modelica_LinearSystems2.Controllers.ADconverter aDconverter(
    y_max=1000,
    y_min=-1000,
    bits=0)
    annotation (Placement(transformation(extent={{-20,-190},{0,-170}})));
  Modelica_LinearSystems2.Controllers.DAconverter dAconverter(
    y_max=1000,
    y_min=-1000,
    bits=0)
    annotation (Placement(transformation(extent={{-20,-220},{0,-200}})));
  Modelica_LinearSystems2.Controllers.Noise noise(y_min=0, y_max=1)
    annotation (Placement(transformation(extent={{-20,-250},{0,-230}})));
equation
  connect(step.y, stateSpace.u[1])      annotation (Line(
      points={{-59,0},{-40,0},{-40,230},{-22,230}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(transferFunction.u, step.y) annotation (Line(
      points={{-22,200},{-40,200},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.u, step.y) annotation (Line(
      points={{-22,170},{-40,170},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter.u, step.y) annotation (Line(
      points={{-22,140},{-40,140},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter1.u, step.y) annotation (Line(
      points={{-22,110},{-40,110},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator.u, step.y) annotation (Line(
      points={{-22,80},{-40,80},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(derivative.u, step.y) annotation (Line(
      points={{-22,50},{-40,50},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder.u, step.y) annotation (Line(
      points={{-22,20},{-40,20},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(secondOrder.u, step.y) annotation (Line(
      points={{-22,-10},{-40,-10},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pI.u, step.y) annotation (Line(
      points={{-22,-40},{-40,-40},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pID.u, step.y) annotation (Line(
      points={{-22,-70},{-40,-70},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID.u_s, step.y) annotation (Line(
      points={{-22,-100},{-40,-100},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.u, step.y) annotation (Line(
      points={{-22,-150},{-40,-150},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(aDconverter.u, step.y) annotation (Line(
      points={{-22,-180},{-40,-180},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dAconverter.u, step.y) annotation (Line(
      points={{-22,-210},{-40,-210},{-40,0},{-59,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(unitDelay.y, PID.u_m) annotation (Line(
      points={{1,-150},{20,-150},{20,-120},{-10,-120},{-10,-112}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (    experiment(StopTime=5));
end TestComponents;
