within Modelica_LinearSystems2.Examples.Utilities;
model DoublePendulum "Double pendulum system"

  parameter Modelica.SIunits.Mass m_trolley = 5 "Mass of trolley";
  parameter Modelica.SIunits.Mass m_load = 20 "Mass of load";
  parameter Modelica.SIunits.Length length = 2 "Length";
  parameter Modelica.SIunits.Angle phi1_start = -80.0/180*pi;
  parameter Modelica.SIunits.Angle phi2_start = 10;
  parameter Modelica.SIunits.AngularVelocity w1_start = 0.0;
  parameter Modelica.SIunits.AngularVelocity w2_start = 0.0;

  constant Real pi = Modelica.Constants.pi;

  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-160,-20},{-120,20}})));
  Modelica.Blocks.Interfaces.RealOutput s
    annotation (Placement(transformation(extent={{120,90},{140,110}})));
  Modelica.Blocks.Interfaces.RealOutput v
    annotation (Placement(transformation(extent={{120,50},{140,70}})));
  Modelica.Blocks.Interfaces.RealOutput phi
    annotation (Placement(transformation(extent={{120,10},{140,30}})));
  Modelica.Blocks.Interfaces.RealOutput w
    annotation (Placement(transformation(extent={{120,-30},{140,-10}})));
  Modelica.Blocks.Interfaces.RealOutput phi1
    annotation (Placement(transformation(extent={{120,-70},{140,-50}})));
  Modelica.Blocks.Interfaces.RealOutput w1
    annotation (Placement(transformation(extent={{120,-110},{140,-90}})));

  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor sensorRelativeTransl(
    get_r_rel=true,
    get_v_rel=true,
    resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameAB.frame_a,
    animation=false) annotation (Placement(transformation(extent={{-90,70},{-70,50}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor sensorRelativeRot1(
    animation=false,
    resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameAB.frame_a,
    get_w_rel=true,
    get_angles=true,
    sequence={1,2,3},
    guessAngle1=0) annotation (Placement(transformation(extent={{-30,70},{-10,50}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor sensorRelativeRot2(
    animation=false,
    resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameAB.frame_a,
    get_w_rel=true,
    get_angles=true,
    sequence={1,2,3},
    guessAngle1=0) annotation (Placement(transformation(extent={{30,-20},{50,0}})));
  inner Modelica.Mechanics.MultiBody.World world(gravityType=Modelica.Mechanics.MultiBody.Types.GravityTypes.
        UniformGravity, animateWorld=false)
                        annotation (Placement(transformation(extent={{-120,10},{-100,30}},
                      rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic jointPrismatic(useAxisFlange=true) annotation (Placement(transformation(extent={{-90,30},{-70,10}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute jointRevolute1(
    n={0,0,1},
    useAxisFlange=true,
    phi(fixed=true, start=phi1_start),
    w(fixed=true, start=w1_start)) annotation (Placement(transformation(extent={{-30,30},{-10,10}},rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute jointRevolute2(
    phi(fixed=true, start=phi2_start),
    w(fixed=true, start=w2_start),
    cylinderDiameter=3*world.defaultJointWidth,
    cylinderColor={0,0,200}) annotation (Placement(transformation(extent={{30,10},{50,30}},rotation=0)));
  Modelica.Mechanics.Translational.Components.Damper damperTranslational(d=0) annotation (Placement(transformation(extent={{-90,-20},{-70,0}})));
  Modelica.Mechanics.Rotational.Components.Damper damperRotational(d=0) annotation (Placement(transformation(extent={{-30,-22},{-10,-2}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyShape bodyTrolley(
    shapeType="box",
    animateSphere=true,
    m=m_trolley,
    sphereDiameter=world.defaultBodyDiameter,
    r={0,0,0},
    r_CM={0,0,0}) annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyPendulum1(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900) annotation (Placement(transformation(extent={{0,10},{20,30}})));
  Modelica.Mechanics.MultiBody.Parts.Body bodyPendulum2(
    m=m_load,
    r_CM={0,0,0},
    specularCoefficient=4*world.defaultSpecularCoefficient,
    sphereDiameter=1.5*world.defaultBodyDiameter) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={70,-30})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyPendulum2Cylinder(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={70,0})));
  Modelica.Mechanics.Translational.Sources.Force force
    annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
protected
  Modelica.Blocks.Sources.Constant const1(k=0.5*Modelica.Constants.pi) annotation (Placement(transformation(extent={{20,54},{32,66}})));
  Modelica.Blocks.Sources.Constant const2(k=0)
    annotation (Placement(transformation(extent={{54,-80},{66,-68}})));
  Modelica.Blocks.Math.Add add1 annotation (Placement(transformation(extent={{50,50},{70,70}})));
  Modelica.Blocks.Math.Add add2
    annotation (Placement(transformation(extent={{80,-70},{100,-50}})));
equation
  connect(damperRotational.flange_b, jointRevolute1.axis) annotation (Line(points={{-10,-12},{-10,10},{-20,10}}, color={0,0,0}));
  connect(jointRevolute1.support, damperRotational.flange_a) annotation (Line(points={{-26,10},{-26,10},{-30,10},{-30,-12}}, color={0,0,0}));
  connect(bodyTrolley.frame_b, jointRevolute1.frame_a) annotation (Line(
      points={{-40,20},{-30,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(jointPrismatic.frame_a, world.frame_b) annotation (Line(
      points={{-90,20},{-100,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(force.flange, jointPrismatic.axis) annotation (Line(
      points={{-70,-40},{-70,14},{-72,14}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damperTranslational.flange_a, jointPrismatic.support) annotation (Line(
      points={{-90,-10},{-90,14},{-84,14}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damperTranslational.flange_b, jointPrismatic.axis) annotation (Line(
      points={{-70,-10},{-70,14},{-72,14}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(jointPrismatic.frame_b, bodyTrolley.frame_a) annotation (Line(
      points={{-70,20},{-60,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(sensorRelativeRot1.frame_b, jointRevolute1.frame_b) annotation (Line(
      points={{-10,60},{-10,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(sensorRelativeRot1.frame_a, jointRevolute1.frame_a) annotation (Line(
      points={{-30,60},{-30,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(u, force.f) annotation (Line(
      points={{-140,0},{-110,0},{-110,-40},{-92,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sensorRelativeRot1.w_rel[3], w) annotation (Line(
      points={{-14,71.6667},{-14,74},{0,74},{0,40},{90,40},{90,-20},{130,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sensorRelativeTransl.v_rel[1], v) annotation (Line(
      points={{-86,70.3333},{-86,90},{110,90},{110,60},{130,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sensorRelativeTransl.r_rel[1], s) annotation (Line(
      points={{-90,70.3333},{-90,100},{130,100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(phi, phi) annotation (Line(
      points={{130,20},{130,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.y, phi) annotation (Line(
      points={{71,60},{100,60},{100,20},{130,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u1, sensorRelativeRot1.angles[3]) annotation (Line(
      points={{48,66},{40,66},{40,78},{-18,78},{-18,71.6667}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sensorRelativeRot2.frame_a, jointRevolute2.frame_a) annotation (Line(
      points={{30,-10},{30,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(sensorRelativeRot2.frame_b, jointRevolute2.frame_b) annotation (Line(
      points={{50,-10},{50,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(const2.y,add2. u2) annotation (Line(
      points={{66.6,-74},{74,-74},{74,-66},{78,-66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add2.u1, sensorRelativeRot2.angles[3]) annotation (Line(
      points={{78,-54},{42,-54},{42,-21.6667}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add2.y, phi1) annotation (Line(
      points={{101,-60},{130,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sensorRelativeRot2.w_rel[3], w1) annotation (Line(
      points={{46,-21.6667},{46,-100},{130,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(bodyPendulum2Cylinder.frame_b, bodyPendulum2.frame_a) annotation (Line(
      points={{70,-10},{70,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyPendulum2Cylinder.frame_a, jointRevolute2.frame_b) annotation (Line(
      points={{70,10},{70,20},{50,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyPendulum1.frame_b, jointRevolute2.frame_a) annotation (Line(
      points={{20,20},{30,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyPendulum1.frame_a, jointRevolute1.frame_b) annotation (Line(
      points={{0,20},{-10,20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(sensorRelativeTransl.frame_a, jointPrismatic.frame_a) annotation (Line(
      points={{-90,60},{-90,20}},
      color={95,95,95},
      thickness=0.5));
  connect(sensorRelativeTransl.frame_b, jointPrismatic.frame_b) annotation (Line(
      points={{-70,60},{-70,20}},
      color={95,95,95},
      thickness=0.5));
  connect(const1.y, add1.u2) annotation (Line(points={{32.6,60},{40,60},{40,54},{48,54}}, color={0,0,127}));
  annotation (
    experiment(
      StartTime=1,
      StopTime=10,
      __Dymola_Algorithm="Dassl"),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-120,-100},{120,100}},
        grid={2,2})),
    Documentation(info="<html>
<p>
The physical model is used in Modelica_LinearSystems2.Examples.StateSpace.doublePendulumController 
where it is being linearized an used as a base for linear controller design. 
The results are used to control the crane system in 
Modelica_LinearSystems2.Controller.Examples.DoublePendulum.mo.
</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-120,-100},{120,100}}),
                    graphics={
        Rectangle(
          extent={{-120,100},{120,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-82,22},{82,18}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Rectangle(extent={{-44,54},{0,28}}, lineColor={0,0,0}),
        Ellipse(
          extent={{-40,34},{-28,22}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255},
          lineThickness=0.5),
        Ellipse(
          extent={{-16,34},{-4,22}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Line(
          points={{-18,-16},{10,-62}},
          color={0,0,0},
          smooth=Smooth.None),
        Ellipse(
          extent={{4,-56},{20,-72}},
          lineColor={0,0,0},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-25,44},{-19,38}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{28,46},{4,46}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{34,40},{10,40}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-22,40},{-18,-16}},
          color={0,0,0},
          smooth=Smooth.None),
        Ellipse(
          extent={{-20,-15},{-14,-21}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid)}));
end DoublePendulum;
