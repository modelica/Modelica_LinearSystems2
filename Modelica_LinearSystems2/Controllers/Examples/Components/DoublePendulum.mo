within Modelica_LinearSystems2.Controllers.Examples.Components;
model DoublePendulum "crane trolley system"

  parameter Modelica.Units.SI.Mass m_trolley = 5;
  parameter Modelica.Units.SI.Mass m_load = 20;
  parameter Modelica.Units.SI.Length length = 2;
  parameter Modelica.Units.SI.Angle phi1_start = -80.0/180*pi;
  parameter Modelica.Units.SI.Angle phi2_start = 10;
  parameter Modelica.Units.SI.AngularVelocity w1_start = 0.0;
  parameter Modelica.Units.SI.AngularVelocity w2_start = 0.0;

  constant Real pi = Modelica.Constants.pi;

  inner Modelica.Mechanics.MultiBody.World world(animateWorld=false,
      animateGravity=false)
    annotation (Placement(transformation(extent={{-140,-80},{-120,-60}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
    annotation (Placement(transformation(extent={{-96,0},{-76,20}})));
  Modelica.Mechanics.Translational.Components.Damper damper1(d=0)
    annotation (Placement(transformation(extent={{-96,14},{-76,34}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute rev(n={0,0,1},useAxisFlange=true,
    phi(fixed=true, start=phi1_start),
    w(fixed=true, start=w1_start))
    annotation (Placement(transformation(extent={{-30,0},{-10,20}}, rotation=0)));
  Modelica.Mechanics.Rotational.Components.Damper damper(d=0)
    annotation (Placement(transformation(extent={{-22,40},{-2,60}},rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.Body body(
    m=m_load,
    r_CM={0,0,0},
    specularCoefficient=4*world.defaultSpecularCoefficient,
    sphereDiameter=1.5*world.defaultBodyDiameter)
    annotation (Placement(transformation(extent={{78,0},{98,20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
    shapeType="box",
    m=m_trolley,
    sphereDiameter=world.defaultBodyDiameter,
    r={0,0,0},
    r_CM={0,0,0})
    annotation (Placement(transformation(extent={{-58,-2},{-38,18}})));
  Modelica.Mechanics.Translational.Sources.Force force
    annotation (Placement(transformation(extent={{-98,34},{-78,54}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles
    annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeVelocity relativeVelocity
    annotation (Placement(transformation(extent={{-96,-30},{-76,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition
    annotation (Placement(transformation(extent={{-96,-60},{-76,-40}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-190,-20},{-150,20}})));
  Modelica.Blocks.Interfaces.RealOutput s
    annotation (Placement(transformation(extent={{150,90},{170,110}})));
  Modelica.Blocks.Interfaces.RealOutput v
    annotation (Placement(transformation(extent={{150,50},{170,70}})));
 Modelica.Blocks.Interfaces.RealOutput phi
    annotation (Placement(transformation(extent={{150,10},{170,30}})));
  Modelica.Blocks.Interfaces.RealOutput w
    annotation (Placement(transformation(extent={{150,-30},{170,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
    relativeAngularVelocity
    annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

  Modelica.Blocks.Sources.Constant const(k=0.5*Modelica.Constants.pi)
    annotation (Placement(transformation(extent={{94,-22},{106,-10}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{116,-10},{136,10}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
    phi(fixed=true, start=phi2_start),
    w(fixed=true, start=w2_start),
    cylinderDiameter=3*world.defaultJointWidth,
    cylinderColor={0,0,200}) annotation (Placement(transformation(extent={{24,0},{
            44,20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles1
    annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
    relativeAngularVelocity1
    annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
 Modelica.Blocks.Interfaces.RealOutput phi1
    annotation (Placement(transformation(extent={{150,-70},{170,-50}})));
  Modelica.Blocks.Interfaces.RealOutput w1
    annotation (Placement(transformation(extent={{150,-110},{170,-90}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{88,-50},{108,-30}})));
  Modelica.Blocks.Sources.Constant const1(k=0)
    annotation (Placement(transformation(extent={{66,-62},{78,-50}})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900)
    annotation (Placement(transformation(extent={{-4,0},{16,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900)
    annotation (Placement(transformation(extent={{52,0},{72,20}})));
equation
  connect(damper.flange_b, rev.axis) annotation (Line(points={{-2,50},{0,50},{0,
          24},{0,20},{-20,20}}, color={0,0,0}));
  connect(rev.support, damper.flange_a) annotation (Line(points={{-26,20},{-26,
          26},{-36,26},{-36,50},{-22,50}}, color={0,0,0}));
  connect(bodyShape.frame_b, rev.frame_a) annotation (Line(
      points={{-38,8},{-34,8},{-34,10},{-30,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(prismatic.frame_a, world.frame_b) annotation (Line(
      points={{-96,10},{-110,10},{-110,-70},{-120,-70}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(force.flange, prismatic.axis) annotation (Line(
      points={{-78,44},{-78,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damper1.flange_a, prismatic.support) annotation (Line(
      points={{-96,24},{-96,16},{-90,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damper1.flange_b, prismatic.axis) annotation (Line(
      points={{-76,24},{-78,24},{-78,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
      points={{-76,10},{-68,10},{-68,8},{-58,8}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeVelocity.frame_b, prismatic.frame_b) annotation (Line(
      points={{-76,-20},{-76,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeVelocity.frame_a, prismatic.frame_a) annotation (Line(
      points={{-96,-20},{-96,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativePosition.frame_b, relativeVelocity.frame_b) annotation (Line(
      points={{-76,-50},{-76,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativePosition.frame_a, relativeVelocity.frame_a) annotation (Line(
      points={{-96,-50},{-96,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles.frame_b, rev.frame_b) annotation (Line(
      points={{-10,-20},{-10,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles.frame_a, rev.frame_a) annotation (Line(
      points={{-30,-20},{-30,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(u, force.f) annotation (Line(
      points={{-170,0},{-136,0},{-136,44},{-100,44}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeAngularVelocity.frame_a, relativeAngles.frame_a) annotation (
      Line(
      points={{-30,-50},{-30,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngularVelocity.frame_b, relativeAngles.frame_b) annotation (
      Line(
      points={{-10,-50},{-10,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngularVelocity.w_rel[3], w) annotation (Line(
      points={{-20,-60.3333},{-20,-66},{120,-66},{120,-20},{160,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeVelocity.v_rel[1], v) annotation (Line(
      points={{-86,-31.6667},{-104,-31.6667},{-104,-32},{-118,-32},{-118,62},{
          42,62},{42,60},{160,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativePosition.r_rel[1], s) annotation (Line(
      points={{-86,-61.6667},{-104,-61.6667},{-104,-58},{-120,-58},{-120,100},{
          160,100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(phi, phi) annotation (Line(
      points={{160,20},{160,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, phi) annotation (Line(
      points={{137,0},{148,0},{148,20},{160,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, add.u2) annotation (Line(
      points={{106.6,-16},{110,-16},{110,-6},{114,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, relativeAngles.angles[3]) annotation (Line(
      points={{114,6},{108,6},{108,-4},{58,-4},{58,-36},{-20,-36},{-20,-30.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeAngles1.frame_a, revolute2.frame_a) annotation (Line(
      points={{24,-20},{24,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles1.frame_b, revolute2.frame_b) annotation (Line(
      points={{44,-20},{44,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles1.frame_a, relativeAngularVelocity1.frame_a)
    annotation (Line(
      points={{24,-20},{24,-50}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngularVelocity1.frame_b, relativeAngles1.frame_b)
    annotation (Line(
      points={{44,-50},{44,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(const1.y, add1.u2)
    annotation (Line(
      points={{78.6,-56},{82,-56},{82,-46},{86,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u1, relativeAngles1.angles[3]) annotation (Line(
      points={{86,-34},{60,-34},{60,-30.3333},{34,-30.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.y, phi1) annotation (Line(
      points={{109,-40},{136,-40},{136,-60},{160,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeAngularVelocity1.w_rel[3], w1) annotation (Line(
      points={{34,-60.3333},{36,-60.3333},{36,-100},{160,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
      points={{72,10},{78,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
      points={{52,10},{44,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
      points={{16,10},{24,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder.frame_a, rev.frame_b) annotation (Line(
      points={{-4,10},{-10,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-150,-100},{150,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
 
Model of a simple double pendulum system. <br>
The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.doublePendulumController where it is being
linearized an used as a base for linear controller design. The results are used to control the crane system
in Modelica_LinearSystems2.Controllers.Examples.DoublePendulum.mo
 
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{150,
            100}}), graphics={
        Rectangle(
          extent={{-150,122},{150,-120}},
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
