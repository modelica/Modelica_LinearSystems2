within Modelica_LinearSystems2.Controller.Examples.Components;
model InverseDoublePendulum "Inverse double pendulum"

  parameter Modelica.SIunits.Mass m_trolley = 1 "Mass of trolley";
  parameter Modelica.SIunits.Mass m_load = 1 "Mass of load on 2nd arm";
  parameter Modelica.SIunits.Length length = 1
    "Total length of double pendulum (i.e. length of each arm = length/2)";
  parameter Modelica.SIunits.Angle phi1_start=90.0/180*pi
    "Initial rotation angle of 1st arm relative to trolley";
  parameter Modelica.SIunits.Angle phi2_start = 0
    "Initial rotation angle of 2nd arm relative to 1st arm";
  parameter Modelica.SIunits.AngularVelocity w1_start = 0.0
    "Initial angular velocity of 1st arm relative to trolley";
  parameter Modelica.SIunits.AngularVelocity w2_start = 0.0
    "Initial angular velocity of 2nd arm relative to 1st arm";

  parameter Modelica.SIunits.Position s_start = 0.0
    "Initial position of trolley relative to world";
  parameter Modelica.SIunits.Velocity v_start = 0.0
    "Initial velocity of trolley relative to world";

  parameter Boolean cartDisturbance=false
    "True, if cart disturbance should be enabled";
  parameter Boolean bodyDisturbance=false
    "True, if body disturbance should be enabled";

  constant Real pi = Modelica.Constants.pi;

  inner Modelica.Mechanics.MultiBody.World world(gravityType=Modelica.Mechanics.MultiBody.Types.GravityTypes.
        UniformGravity,
    animateWorld=false,
    animateGravity=false)
                        annotation (Placement(transformation(extent={{-140,-80},
            {-120,-60}},
                      rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true,
    s(start=s_start, fixed=true),
    animation=false,
    v(start=v_start))
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica.Mechanics.Translational.Components.Damper damper1(d=0)
    annotation (Placement(transformation(extent={{-100,18},{-80,38}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(
    n={0,0,1},
    useAxisFlange=true,
    phi(fixed=true, start=phi1_start),
    w(fixed=true, start=w1_start),
    cylinderColor=bodyShape.color,
    cylinderDiameter=2*bodyCylinder.diameter)
                               annotation (Placement(transformation(extent={{-30,0},
            {-10,20}},      rotation=0)));
  Modelica.Mechanics.Rotational.Components.Damper damper(d=0)
    annotation (Placement(transformation(extent={{-30,22},{-10,42}},
                                                                   rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.Body body(
    m=m_load,
    r_CM={0,0,0},
    specularCoefficient=4*world.defaultSpecularCoefficient,
    sphereDiameter=1.5*world.defaultBodyDiameter,
    sphereColor=bodyCylinder1.color)
    annotation (Placement(transformation(extent={{80,0},{100,20}},rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
    m=m_trolley,
    sphereDiameter=world.defaultBodyDiameter,
    animateSphere=false,
    shapeType="box",
    lengthDirection={0,-1,0},
    widthDirection={1,0,0},
    length=0.1,
    r_shape=0.3*revolute1.cylinderDiameter*{0,-1,0},
    width=0.3,
    height=0.5*bodyShape.width,
    color={0,0,0})
    annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
  Modelica.Mechanics.Translational.Sources.Force force
    annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles
    annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeVelocity relativeVelocity
    annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition
    annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
    relativeAngularVelocity
    annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

  Modelica.Blocks.Sources.Constant const(k=-0.5*Modelica.Constants.pi)
    annotation (Placement(transformation(extent={{98,34},{110,46}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{120,30},{140,10}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
    phi(fixed=true, start=phi2_start),
    w(fixed=true, start=w2_start),
    specularCoefficient=0.7,
    cylinderDiameter=2*bodyCylinder.diameter,
    cylinderColor=bodyCylinder.color)                    annotation (Placement(transformation(extent={{24,0},{
            44,20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles1
    annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
  Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
    relativeAngularVelocity1
    annotation (Placement(transformation(extent={{24,-64},{44,-44}})));
  Modelica.Blocks.Sources.Constant const1(k=0)
    annotation (Placement(transformation(extent={{98,-82},{110,-70}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{120,-70},{140,-50}})));
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
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
    r={length/2,0,0},
    specularCoefficient=0.7,
    diameter=0.03,
    density=1000,
    color={0,128,0})
    annotation (Placement(transformation(extent={{-2,0},{18,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
    r={length/2,0,0},
    specularCoefficient=0.7,
    diameter=0.03,
    density=1000,
    color={0,0,255})
    annotation (Placement(transformation(extent={{52,0},{72,20}})));
  Modelica.Blocks.Interfaces.RealInput dist if cartDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-80,120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-60,122})));
  Modelica.Mechanics.Translational.Sources.Force distrubanceForceCart if cartDisturbance
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-80,80})));
  Modelica.Blocks.Interfaces.RealInput dist2 if bodyDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={80,120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={60,122})));
  Modelica.Mechanics.MultiBody.Forces.Torque torque if bodyDisturbance
    annotation (Placement(transformation(extent={{40,60},{60,80}})));
  Modelica.Blocks.Sources.Constant const2[2](k={0,0}) if bodyDisturbance
    annotation (Placement(transformation(extent={{0,70},{20,90}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape1(
    shapeType="cylinder",
    lengthDirection=revolute2.n,
    widthDirection={0,0,1},
    length=1.2*revolute2.cylinderLength,
    width=bodyCylinder1.diameter,
    height=bodyCylinder1.diameter,
    color=bodyCylinder1.color,
    r_shape=-0.5*fixedShape1.length*fixedShape1.lengthDirection)
    annotation (Placement(transformation(extent={{52,20},{72,40}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
    shapeType="cylinder",
    lengthDirection=revolute1.n,
    widthDirection={0,0,1},
    length=1.2*revolute1.cylinderLength,
    width=bodyCylinder.diameter,
    height=bodyCylinder.diameter,
    color=bodyCylinder.color,
    r_shape=-0.5*fixedShape.length*fixedShape.lengthDirection)
    annotation (Placement(transformation(extent={{0,20},{20,40}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape2(
    lengthDirection={0,-1,0},
    widthDirection={1,0,0},
    height=0.5*bodyShape.height,
    color={100,100,100},
    length=0.5*bodyShape.length,
    width=50*bodyShape.width,
    r_shape=bodyShape.r_shape + 0.5*(bodyShape.length - fixedShape2.length)*{200,
        -1,0})
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
equation
  connect(damper.flange_b, revolute1.axis)
                                     annotation (Line(points={{-10,32},{-8,32},
          {-8,22},{-8,20},{-20,20}},
                                  color={0,0,0}));
  connect(revolute1.support, damper.flange_a)
                                        annotation (Line(points={{-26,20},{-26,
          20},{-34,20},{-34,32},{-30,32}}, color={0,0,0}));
  connect(bodyShape.frame_b, revolute1.frame_a)
                                          annotation (Line(
      points={{-38,10},{-30,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(prismatic.frame_a, world.frame_b) annotation (Line(
      points={{-100,10},{-110,10},{-110,-70},{-120,-70}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(force.flange, prismatic.axis) annotation (Line(
      points={{-80,50},{-70,50},{-70,16},{-82,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damper1.flange_a, prismatic.support) annotation (Line(
      points={{-100,28},{-100,16},{-94,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(damper1.flange_b, prismatic.axis) annotation (Line(
      points={{-80,28},{-80,16},{-82,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
      points={{-80,10},{-58,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeVelocity.frame_b, prismatic.frame_b) annotation (Line(
      points={{-80,-20},{-80,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeVelocity.frame_a, prismatic.frame_a) annotation (Line(
      points={{-100,-20},{-100,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativePosition.frame_b, relativeVelocity.frame_b) annotation (Line(
      points={{-80,-50},{-80,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativePosition.frame_a, relativeVelocity.frame_a) annotation (Line(
      points={{-100,-50},{-100,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles.frame_b, revolute1.frame_b)
                                               annotation (Line(
      points={{-10,-20},{-10,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngles.frame_a, revolute1.frame_a)
                                               annotation (Line(
      points={{-30,-20},{-30,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(u, force.f) annotation (Line(
      points={{-170,0},{-136,0},{-136,50},{-102,50}},
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
      points={{-20,-60.3333},{-20,-72},{80,-72},{80,-30},{140,-30},{140,-20},{
          160,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeVelocity.v_rel[1], v) annotation (Line(
      points={{-90,-31.6667},{-108,-31.6667},{-108,-34},{-122,-34},{-122,60},{
          160,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativePosition.r_rel[1], s) annotation (Line(
      points={{-90,-61.6667},{-108,-61.6667},{-108,-58},{-124,-58},{-124,100},{
          160,100}},
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
  connect(add.u1, relativeAngles.angles[3]) annotation (Line(
      points={{118,14},{114,14},{114,-20},{60,-20},{60,-36},{-20,-36},{-20,
          -30.3333}},
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
      points={{24,-20},{24,-54}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativeAngularVelocity1.frame_b, relativeAngles1.frame_b)
    annotation (Line(
      points={{44,-54},{44,-20}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(const1.y, add1.u2)
                           annotation (Line(
      points={{110.6,-76},{114,-76},{114,-66},{118,-66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.u1, relativeAngles1.angles[3]) annotation (Line(
      points={{118,-54},{70,-54},{70,-34},{34,-34},{34,-30.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add1.y, phi1) annotation (Line(
      points={{141,-60},{160,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativeAngularVelocity1.w_rel[3], w1) annotation (Line(
      points={{34,-64.3333},{34,-100},{160,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
      points={{18,10},{24,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder.frame_a, revolute1.frame_b)
                                             annotation (Line(
      points={{-2,10},{-10,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
      points={{72,10},{80,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
      points={{52,10},{44,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(distrubanceForceCart.f, dist) annotation (Line(
      points={{-80,92},{-80,120}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(distrubanceForceCart.flange, prismatic.axis) annotation (Line(
      points={{-80,70},{-70,70},{-70,16},{-82,16}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(dist2, torque.torque[3]) annotation (Line(
      points={{80,120},{80,90},{44,90},{44,83.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const2.y, torque.torque[1:2]) annotation (Line(
      points={{21,80},{34,80},{34,86},{44,86},{44,80.6667}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(torque.frame_a, revolute2.frame_a) annotation (Line(
      points={{40,70},{40,70},{24,70},{24,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(torque.frame_b, bodyCylinder1.frame_b) annotation (Line(
      points={{60,70},{74,70},{74,10},{72,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape1.frame_a, revolute2.frame_b) annotation (Line(
      points={{52,30},{48,30},{48,10},{44,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape.frame_a, revolute1.frame_b) annotation (Line(
      points={{0,30},{-6,30},{-6,10},{-10,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape2.frame_a, world.frame_b) annotation (Line(
      points={{-100,-90},{-110,-90},{-110,-70},{-120,-70}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=50, Algorithm="Dassl"),
    Documentation(info="<html>

Model of a simple inverted double pendulum system. The mdel is the same as in Modelica_Controller.Examples.Components.DoublePendulum but with different initial values because the initial values are used as a working point for linearization.<br>
The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.inverseDoublePendulumController where it is being
linearized an used as a base for linear controller design. The results are used to control the crane system
in Modelica_Controller.Examples.InverseDoublePendulum.mo

</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,102},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-82,-74},{82,-78}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Rectangle(extent={{-44,-42},{0,-68}}, lineColor={0,0,0}),
        Ellipse(
          extent={{-40,-62},{-28,-74}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255},
          lineThickness=0.5),
        Ellipse(
          extent={{-16,-62},{-4,-74}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Line(
          points={{-16,-6},{-22,-56}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Ellipse(
          extent={{-32,54},{-20,42}},
          lineColor={0,0,0},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-25,-52},{-19,-58}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{34,-56},{10,-56}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{28,-64},{4,-64}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-34,54},{-34,54},{-36,52},{-36,46},{-34,44}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-38,56},{-38,56},{-40,54},{-40,48},{-38,46}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-26,42},{-16,-6}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Ellipse(
          extent={{-20,-2},{-14,-8}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}));
end InverseDoublePendulum;
