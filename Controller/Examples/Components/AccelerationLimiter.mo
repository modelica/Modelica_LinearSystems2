within Modelica_LinearSystems2.Controller.Examples.Components;
block AccelerationLimiter
  "output follows input with limited acceleration and optionally limited velocity"
 extends Interfaces.PartialSampledBlock;
  parameter Boolean velocityLimitation = true "true if velocity to be limited"   annotation(Evaluate=true,choices(checkBox=true));
  parameter Real v_limit(min=Modelica.Constants.eps)=1
    "maximal absolute velocity" annotation(Dialog(enable=velocityLimitation));
  parameter Real a_limit(min=Modelica.Constants.eps)=1
    "maximal absolute acceleration";

  parameter Real y1_start=0 "start value of integrtor1" annotation(Evaluate=true,Dialog(tab="Advanced options",group = "integrator1"));
  parameter Boolean withDelay1=blockType==Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Discrete
    "delays the input of integrator1" annotation(Evaluate=true,Dialog(tab="Advanced options",group = "integrator1",
                enable=blockType==Modelica_Controller.Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Discrete));
  parameter Real y2_start=0 "start value of integrtor2" annotation(Evaluate=true,Dialog(tab="Advanced options",group = "integrator2"));
  parameter Boolean withDelay2=blockType==Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Discrete
    "delays the input of integrator2"  annotation(Dialog(tab="Advanced options",group = "integrator2",
                enable=blockType==Modelica_Controller.Types.BlockTypeWithGlobalDefault.Discrete or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Discrete));

  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-85,5},{-75,-5}})));
  Modelica.Blocks.Nonlinear.Limiter limiter1(uMax=v_limit) if velocityLimitation
    annotation (Placement(transformation(extent={{-44,5},{-34,15}})));
  Modelica.Blocks.Math.Gain gain(k=1) if not velocityLimitation
    annotation (Placement(transformation(extent={{-44,-16},{-34,-6}})));
  Integrator integrator1(
    k=1,
    withDelay=withDelay1,
    y_start=y1_start)
    annotation (Placement(transformation(extent={{34,-6},{46,6}})));
  Modelica.Blocks.Math.Feedback feedback1
    annotation (Placement(transformation(extent={{-25,5},{-15,-5}})));
  Modelica.Blocks.Nonlinear.Limiter limiter2(uMax=a_limit)
    annotation (Placement(transformation(extent={{-6,-5},{4,5}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput s
    annotation (Placement(transformation(extent={{100,50},{120,70}})));
  Modelica.Blocks.Interfaces.RealOutput v
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealOutput a
    annotation (Placement(transformation(extent={{100,-70},{120,-50}})));
  Modelica.Blocks.Math.Abs abs1
    annotation (Placement(transformation(extent={{25,35},{15,45}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{5,55},{-5,65}})));
  Modelica.Blocks.Math.Gain gain1(k=1/(2*a_limit))
    annotation (Placement(transformation(extent={{-16,55},{-26,65}})));
  Internal.Add add
    annotation (Placement(transformation(extent={{-55,55},{-65,65}})));
  Integrator integrator2(
    k=1,
    withDelay=withDelay2,
    y_start=y2_start)
    annotation (Placement(transformation(extent={{74,-6},{86,6}})));
equation
  connect(feedback.y, limiter1.u) annotation (Line(
      points={{-75.5,0},{-50,0},{-50,10},{-45,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback1.y, limiter2.u) annotation (Line(
      points={{-15.5,0},{-7,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter2.y, integrator1.u) annotation (Line(
      points={{4.5,0},{32.8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, feedback1.u2) annotation (Line(
      points={{46.6,0},{60,0},{60,20},{-20,20},{-20,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter1.y, feedback1.u1) annotation (Line(
      points={{-33.5,10},{-28,10},{-28,0},{-24,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, feedback.y) annotation (Line(
      points={{-45,-11},{-50,-11},{-50,0},{-75.5,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, feedback1.u1) annotation (Line(
      points={{-33.5,-11},{-28,-11},{-28,0},{-24,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, abs1.u) annotation (Line(
      points={{46.6,0},{60,0},{60,40},{26,40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(abs1.y, product.u2) annotation (Line(
      points={{14.5,40},{10,40},{10,57},{6,57}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product.u1, integrator1.y) annotation (Line(
      points={{6,63},{60,63},{60,0},{46.6,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain1.u, product.y) annotation (Line(
      points={{-15,60},{-5.5,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, gain1.y) annotation (Line(
      points={{-56,60},{-26.5,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, feedback.u2) annotation (Line(
      points={{-64.5,60},{-80,60},{-80,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator1.y, v) annotation (Line(
      points={{46.6,0},{60,0},{60,-20},{96,-20},{96,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(limiter2.y, a) annotation (Line(
      points={{4.5,0},{20,0},{20,-60},{110,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, u) annotation (Line(
      points={{-84,0},{-96,0},{-96,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.u, integrator1.y) annotation (Line(
      points={{72.8,0},{46.6,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.y, add.u1) annotation (Line(
      points={{86.6,0},{90,0},{90,80},{-60,80},{-60,64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(integrator2.y, s) annotation (Line(
      points={{86.6,0},{90,0},{90,60},{110,60}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics), Icon(graphics={
        Text(
          extent={{-38,-36},{92,-80}},
          lineColor={0,0,255},
          textString="acceleration"),
        Text(
          visible=velocityLimitation,
          extent={{8,18},{90,-10}},
          lineColor={0,0,255},
          textString="velocity"),
        Text(
          extent={{-88,86},{82,26}},
          lineColor={0,0,255},
          textString="input with limited"),
        Text(
          visible=velocityLimitation,
          extent={{-106,-14},{-28,-38}},
          lineColor={0,0,255},
          textString="and")}));
end AccelerationLimiter;
