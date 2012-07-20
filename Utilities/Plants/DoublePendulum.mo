within Modelica_LinearSystems2.Utilities.Plants;
model DoublePendulum "Multibody model of crane trolley"

  parameter Modelica.SIunits.Mass m_trolley = 5 "Mass of trolley";
  parameter Modelica.SIunits.Mass m_load = 20 "Mass of load on 2nd arm";
  parameter Modelica.SIunits.Length length = 2
    "Total length of double pendulum (i.e. length of each arm = length/2)";

  parameter Modelica.SIunits.Position s_start = 0.0
    "Initial position of trolley relative to world";
  parameter Modelica.SIunits.Velocity v_start = 0.0
    "Initial velocity of trolley relative to world";
  parameter Modelica.SIunits.Angle phi1_start = -80.0/180*pi
    "Initial rotation angle of 1st arm relative to trolley";
  parameter Modelica.SIunits.Angle phi2_start = 10
    "Initial rotation angle of 2nd arm relative to 1st arm";
  parameter Modelica.SIunits.AngularVelocity w1_start = 0.0
    "Initial angular velocity of 1st arm relative to trolley";
  parameter Modelica.SIunits.AngularVelocity w2_start = 0.0
    "Initial angular velocity of 2nd arm relative to 1st arm";

  parameter Modelica.SIunits.RotationalDampingConstant d=0
    "Damping constant for revolute joint of 1st arm ";

  constant Real pi = Modelica.Constants.pi;

  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-190,-20},{-150,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput s
    annotation (Placement(transformation(extent={{150,90},{170,110}}),
        iconTransformation(extent={{100,90},{120,110}})));
  Modelica.Blocks.Interfaces.RealOutput v
    annotation (Placement(transformation(extent={{150,50},{170,70}}),
        iconTransformation(extent={{100,50},{120,70}})));
 Modelica.Blocks.Interfaces.RealOutput phi
    annotation (Placement(transformation(extent={{150,10},{170,30}}),
        iconTransformation(extent={{100,10},{120,30}})));
  Modelica.Blocks.Interfaces.RealOutput w
    annotation (Placement(transformation(extent={{150,-30},{170,-10}}),
        iconTransformation(extent={{100,-30},{120,-10}})));

 Modelica.Blocks.Interfaces.RealOutput phi1
    annotation (Placement(transformation(extent={{150,-70},{170,-50}}),
        iconTransformation(extent={{100,-70},{120,-50}})));
  Modelica.Blocks.Interfaces.RealOutput w1
    annotation (Placement(transformation(extent={{150,-110},{170,-90}}),
        iconTransformation(extent={{100,-110},{120,-90}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor relativeSensorPrismatic(
    animation=false,
    get_r_rel=true,
    get_v_rel=true)
    annotation (Placement(transformation(extent={{-100,50},{-80,30}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor relativeSensorAng1(
    get_w_rel=true,
    get_angles=true,
    animation=false)
    annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeSensor relativeSensorAng2(
    get_w_rel=true,
    get_angles=true,
    animation=false)
    annotation (Placement(transformation(extent={{20,-30},{40,-10}})));
  inner Modelica.Mechanics.MultiBody.World world(animateWorld=false,
      animateGravity=false)
                        annotation (Placement(transformation(extent={{-130,0},{-110,
            20}},     rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true,
    s(fixed=true, start=s_start),
    v(fixed=true, start=v_start))
    annotation (Placement(transformation(extent={{-100,20},{-80,0}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(
    n={0,0,1},
    useAxisFlange=true,
    phi(fixed=true, start=phi1_start),
    w(fixed=true, start=w1_start))
                               annotation (Placement(transformation(extent={{-40,0},
            {-20,20}},      rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
    phi(fixed=true, start=phi2_start),
    w(fixed=true, start=w2_start),
    cylinderDiameter=3*world.defaultJointWidth,
    cylinderColor={0,0,200})                             annotation (Placement(transformation(extent={{20,0},{
            40,20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
    shapeType="box",
    m=m_trolley,
    sphereDiameter=world.defaultBodyDiameter,
    r={0,0,0},
    r_CM={0,0,0},
    useQuaternions=false)
    annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
  Modelica.Mechanics.MultiBody.Parts.Body body(
    m=m_load,
    r_CM={0,0,0},
    specularCoefficient=4*world.defaultSpecularCoefficient,
    sphereDiameter=1.5*world.defaultBodyDiameter)
    annotation (Placement(transformation(extent={{80,0},{100,20}},rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900,
    useQuaternions=false)
    annotation (Placement(transformation(extent={{-10,0},{10,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
    r={length/2,0,0},
    specularCoefficient=0.7,
    color={0,0,0},
    diameter=0.05,
    density=900,
    useQuaternions=false)
    annotation (Placement(transformation(extent={{50,0},{70,20}})));
  Modelica.Mechanics.Translational.Components.Damper damperTrans1D(d=0)
    annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
  Modelica.Mechanics.Rotational.Components.Damper damperRot1D(d=d)
    annotation (Placement(transformation(extent={{-40,20},{-20,40}},
                                                                   rotation=0)));
  Modelica.Mechanics.Translational.Sources.Force force
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}})));
  Modelica.Blocks.Sources.Constant const(k=0.5*Modelica.Constants.pi)
    annotation (Placement(transformation(extent={{98,34},{110,46}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{120,30},{140,10}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{120,-70},{140,-50}})));
  Modelica.Blocks.Sources.Constant const1(k=0)
    annotation (Placement(transformation(extent={{98,-82},{110,-70}})));
equation
  connect(damperRot1D.flange_b, revolute1.axis)
                                     annotation (Line(points={{-20,30},{-20,30},
          {-20,22},{-20,20},{-30,20}},
                                  color={0,0,0}));
  connect(revolute1.support, damperRot1D.flange_a)
                                        annotation (Line(points={{-36,20},{-36,20},
          {-40,20},{-40,30}},              color={0,0,0}));
  connect(bodyShape.frame_b, revolute1.frame_a)
                                          annotation (Line(
      points={{-50,10},{-40,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(force.flange, prismatic.axis) annotation (Line(
      points={{-80,-40},{-80,4},{-82,4}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damperTrans1D.flange_a, prismatic.support)
                                               annotation (Line(
      points={{-100,-20},{-100,4},{-94,4}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damperTrans1D.flange_b, prismatic.axis)
                                            annotation (Line(
      points={{-80,-20},{-80,4},{-82,4}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
      points={{-80,10},{-70,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(u, force.f) annotation (Line(
      points={{-170,0},{-140,0},{-140,-40},{-102,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeSensorAng1.w_rel[3], w) annotation (Line(
      points={{-24,-30.3333},{-24,-40},{140,-40},{140,-20},{160,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeSensorPrismatic.v_rel[1], v) annotation (Line(
      points={{-96,51.6667},{-96,60},{160,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeSensorPrismatic.r_rel[1], s) annotation (Line(
      points={{-100,51.6667},{-100,64},{140,64},{140,100},{160,100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, phi) annotation (Line(
      points={{141,20},{160,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, add.u2) annotation (Line(
      points={{110.6,40},{114,40},{114,26},{118,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, relativeSensorAng1.angles[3]) annotation (Line(
      points={{118,14},{112,14},{112,-36},{-28,-36},{-28,-30.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const1.y, add1.u2)
                           annotation (Line(
      points={{110.6,-76},{114,-76},{114,-66},{118,-66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u1, relativeSensorAng2.angles[3]) annotation (Line(
      points={{118,-54},{32,-54},{32,-30.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.y, phi1) annotation (Line(
      points={{141,-60},{160,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeSensorAng2.w_rel[3], w1) annotation (Line(
      points={{36,-30.3333},{36,-100},{160,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
      points={{70,10},{80,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
      points={{50,10},{40,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
      points={{10,10},{20,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder.frame_a, revolute1.frame_b)
                                             annotation (Line(
      points={{-10,10},{-20,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeSensorAng1.frame_a, revolute1.frame_a)
                                                   annotation (Line(
      points={{-40,-20},{-44,-20},{-44,10},{-40,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeSensorAng1.frame_b, revolute1.frame_b)
                                                   annotation (Line(
      points={{-20,-20},{-16,-20},{-16,10},{-20,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(revolute2.frame_a, relativeSensorAng2.frame_a) annotation (Line(
      points={{20,10},{16,10},{16,-20},{20,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(revolute2.frame_b, relativeSensorAng2.frame_b) annotation (Line(
      points={{40,10},{44,10},{44,-20},{40,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(prismatic.frame_b, relativeSensorPrismatic.frame_b) annotation (Line(
      points={{-80,10},{-76,10},{-76,40},{-80,40}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(world.frame_b, prismatic.frame_a) annotation (Line(
      points={{-110,10},{-100,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeSensorPrismatic.frame_a, prismatic.frame_a) annotation (Line(
      points={{-100,40},{-104,40},{-104,10},{-100,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    experiment(
      StartTime=1,
      StopTime=10,
      Algorithm="Dassl"),
    Documentation(info="<html>
<p>Multibody model of a simple double pendulum system. This physical model is used in various models and functions of the library e.g. for linearization or as a base for linear controller design. </p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
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
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,140},{150,100}},
          lineColor={0,0,255},
          textString="%name")}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{150,
            100}})));
end DoublePendulum;
