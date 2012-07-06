within Modelica_LinearSystems2.Utilities.Plants;
model DoublePendulumInverse "Multibody model of inverse crane trolley"
  extends DoublePendulum(
    m_trolley = 1,
    m_load = 1,
    length = 1,
    phi1_start=90.0/180*pi,
    phi2_start = 0,
    w1_start = 0.0,
    w2_start = 0.0,
    s_start = 0.0,
    v_start = 0.0);

  parameter Boolean cartDisturbance=false
    "True, if cart disturbance should be enabled";
  parameter Boolean bodyDisturbance=false
    "True, if body disturbance should be enabled";

  Modelica.Blocks.Interfaces.RealInput dist if cartDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-80,120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-60,122})));
  Modelica.Blocks.Interfaces.RealInput dist2 if bodyDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={80,120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={60,122})));
  Modelica.Mechanics.Translational.Sources.Force distrubanceForceCart if cartDisturbance
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-80,80})));
  Modelica.Mechanics.MultiBody.Forces.Torque torque if bodyDisturbance
    annotation (Placement(transformation(extent={{20,64},{40,84}})));
  Modelica.Blocks.Sources.Constant constZero(k=0) if bodyDisturbance
    annotation (Placement(transformation(extent={{-20,80},{0,100}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
    shapeType="cylinder",
    lengthDirection=revolute1.n,
    widthDirection={0,0,1},
    length=1.2*revolute1.cylinderLength,
    width=bodyCylinder.diameter,
    height=bodyCylinder.diameter,
    color=bodyCylinder.color,
    r_shape=-0.5*fixedShape.length*fixedShape.lengthDirection)
    annotation (Placement(transformation(extent={{-10,30},{10,50}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape1(
    shapeType="cylinder",
    lengthDirection=revolute2.n,
    widthDirection={0,0,1},
    length=1.2*revolute2.cylinderLength,
    width=bodyCylinder1.diameter,
    height=bodyCylinder1.diameter,
    color=bodyCylinder1.color,
    r_shape=-0.5*fixedShape1.length*fixedShape1.lengthDirection)
    annotation (Placement(transformation(extent={{50,30},{70,50}})));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape2(
    lengthDirection={0,-1,0},
    widthDirection={1,0,0},
    height=0.5*bodyShape.height,
    color={100,100,100},
    length=0.5*bodyShape.length,
    width=50*bodyShape.width,
    r_shape=bodyShape.r_shape + 0.5*(bodyShape.length - fixedShape2.length)*{200,
        -1,0})
    annotation (Placement(transformation(extent={{-110,30},{-130,50}})));
equation
  connect(distrubanceForceCart.f,dist)  annotation (Line(
      points={{-80,92},{-80,120}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(distrubanceForceCart.flange, prismatic.axis) annotation (Line(
      points={{-80,70},{-74,70},{-74,0},{-80,0},{-80,4},{-82,4}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(dist2,torque. torque[3]) annotation (Line(
      points={{80,120},{80,90},{24,90},{24,87.3333}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(torque.frame_a, revolute2.frame_a) annotation (Line(
      points={{20,74},{16,74},{16,10},{20,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(constZero.y, torque.torque[1]) annotation (Line(
      points={{1,90},{24.5,90},{24.5,84.6667},{24,84.6667}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(constZero.y, torque.torque[2]) annotation (Line(
      points={{1,90},{24,90},{24,86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(torque.frame_b, revolute2.frame_b) annotation (Line(
      points={{40,74},{40,74},{44,74},{44,10},{40,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape1.frame_a, revolute2.frame_b) annotation (Line(
      points={{50,40},{44,40},{44,10},{40,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape2.frame_a, prismatic.frame_a) annotation (Line(
      points={{-110,40},{-104,40},{-104,10},{-100,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedShape.frame_a, revolute1.frame_b) annotation (Line(
      points={{-10,40},{-16,40},{-16,10},{-20,10}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  annotation ( Documentation(info="<html>
<p>Multibody model of a simple inverted double pendulum system. This physical model is used in various models and functions of the library e.g. for linearization or as a base for linear controller design. The mdel is the same as in Modelica_Controller.Examples.Components.DoublePendulum but with different initial values because the initial values are used as a working point for linearization.</p>
</html>"));
end DoublePendulumInverse;
