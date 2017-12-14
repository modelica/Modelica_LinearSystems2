within Modelica_LinearSystems2.Controllers.Examples.Components;
block AccelerationLimiter
  "Output follows input with limited acceleration and optionally limited velocity"
  extends Interfaces.PartialSampledBlock;
  parameter Boolean velocityLimitation = true "True, if velocity to be limited"
    annotation(Evaluate=true,choices(checkBox=true));
  parameter Real v_limit(min=Modelica.Constants.eps)=1
    "Maximum absolute velocity" annotation(Dialog(enable=velocityLimitation));
  parameter Real a_limit(min=Modelica.Constants.eps)=1
    "Maximum absolute acceleration";

  parameter Real y1_start=0 "Start value of integrator1" annotation(Evaluate=true,Dialog(tab="Advanced options",group = "Integrator1"));
  parameter Boolean withDelay1=blockType==Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
    sampleClock.blockType == Types.BlockType.Discrete
    "Delays the input of integrator1"
    annotation(Evaluate=true,Dialog(tab="Advanced options",group = "Integrator1",
      enable=blockType==Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Discrete or
        blockType == Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
        sampleClock.blockType == Modelica_LinearSystems2.Controllers.Types.BlockType.Discrete));
  parameter Real y2_start=0 "Start value of integrator2" annotation(Evaluate=true,Dialog(tab="Advanced options",group = "Integrator2"));
  parameter Boolean withDelay2=blockType==Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
    sampleClock.blockType == Types.BlockType.Discrete
    "Delays the input of integrator2"
    annotation(Dialog(tab="Advanced options",group = "Integrator2",
      enable=blockType==Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Discrete or
        blockType == Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
        sampleClock.blockType == Modelica_LinearSystems2.Controllers.Types.BlockType.Discrete));

  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-85,5},{-75,-5}})));
  Modelica.Blocks.Nonlinear.Limiter limiter1(uMax=v_limit) if velocityLimitation
    annotation (Placement(transformation(extent={{-60,5},{-50,15}})));
  Modelica.Blocks.Math.Gain gain(k=1) if not velocityLimitation
    annotation (Placement(transformation(extent={{-60,-16},{-50,-6}})));
  Integrator integrator1(
    k=1,
    withDelay=withDelay1,
    y_start=y1_start)
    annotation (Placement(transformation(extent={{14,-6},{26,6}})));
  Modelica.Blocks.Math.Feedback feedback1
    annotation (Placement(transformation(extent={{-37,5},{-27,-5}})));
  Modelica.Blocks.Nonlinear.Limiter limiter2(uMax=a_limit)
    annotation (Placement(transformation(extent={{-18,-5},{-8,5}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput s
    annotation (Placement(transformation(extent={{100,50},{120,70}})));
  Modelica.Blocks.Interfaces.RealOutput v
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealOutput a
    annotation (Placement(transformation(extent={{100,-70},{120,-50}})));
  Modelica.Blocks.Math.Abs abs1
    annotation (Placement(transformation(extent={{35,35},{25,45}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{15,55},{5,65}})));
  Modelica.Blocks.Math.Gain gain1(k=1/(2*a_limit))
    annotation (Placement(transformation(extent={{-6,55},{-16,65}})));
  Internal.Add add
    annotation (Placement(transformation(extent={{-55,55},{-65,65}})));
  Integrator integrator2(
    k=1,
    withDelay=withDelay2,
    y_start=y2_start)
    annotation (Placement(transformation(extent={{54,-6},{66,6}})));
equation
  connect(feedback.y, limiter1.u) annotation (Line(
      points={{-75.5,0},{-66,0},{-66,10},{-61,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback1.y, limiter2.u) annotation (Line(
      points={{-27.5,0},{-19,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter2.y, integrator1.u) annotation (Line(
      points={{-7.5,0},{12.8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, feedback1.u2) annotation (Line(
      points={{26.6,0},{40,0},{40,20},{-32,20},{-32,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter1.y, feedback1.u1) annotation (Line(
      points={{-49.5,10},{-44,10},{-44,0},{-36,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, feedback.y) annotation (Line(
      points={{-61,-11},{-66,-11},{-66,0},{-75.5,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, feedback1.u1) annotation (Line(
      points={{-49.5,-11},{-44,-11},{-44,0},{-36,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, abs1.u) annotation (Line(
      points={{26.6,0},{40,0},{40,40},{36,40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(abs1.y, product.u2) annotation (Line(
      points={{24.5,40},{20,40},{20,57},{16,57}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product.u1, integrator1.y) annotation (Line(
      points={{16,63},{40,63},{40,0},{26.6,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain1.u, product.y) annotation (Line(
      points={{-5,60},{4.5,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, gain1.y) annotation (Line(
      points={{-56,60},{-16.5,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, feedback.u2) annotation (Line(
      points={{-64.5,60},{-80,60},{-80,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, v) annotation (Line(
      points={{26.6,0},{40,0},{40,-20},{90,-20},{90,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter2.y, a) annotation (Line(
      points={{-7.5,0},{0,0},{0,-60},{110,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, u) annotation (Line(
      points={{-84,0},{-96,0},{-96,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.u, integrator1.y) annotation (Line(
      points={{52.8,0},{26.6,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.y, add.u1) annotation (Line(
      points={{66.6,0},{80,0},{80,80},{-60,80},{-60,64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.y, s) annotation (Line(
      points={{66.6,0},{80,0},{80,60},{110,60}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation ( Icon(graphics={
        Text(
          extent={{-38,-36},{92,-80}},
          lineColor={0,0,0},
          textString="acceleration"),
        Text(
          visible=velocityLimitation,
          extent={{8,18},{90,-10}},
          lineColor={0,0,0},
          textString="velocity"),
        Text(
          extent={{-88,86},{82,26}},
          lineColor={0,0,0},
          textString="input with limited"),
        Text(
          visible=velocityLimitation,
          extent={{-106,-14},{-28,-38}},
          lineColor={0,0,0},
          textString="and")}), Diagram(graphics));
end AccelerationLimiter;
