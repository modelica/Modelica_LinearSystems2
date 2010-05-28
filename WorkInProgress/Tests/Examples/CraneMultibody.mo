within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model CraneMultibody
  import SI = Modelica.SIunits;
  import Modelica.Math.*;
  parameter SI.Mass m_crab=1000 "mass of crab";
  parameter SI.Mass m_load=4000 "mass of load";
  parameter SI.Length l=10 "length of rope";
  parameter SI.Acceleration g = 9.81 "Gravity acceleration";
  parameter Real d=1e6 "damping";
  parameter Real J=100;
  Modelica.Blocks.Interfaces.RealInput force 
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
      Modelica.Blocks.Interfaces.RealOutput x1 "Angle of pendulum" 
    annotation (Placement(transformation(extent={{100,42},{120,62}},
          rotation=0), iconTransformation(extent={{100,80},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput x2 "Angular velocity of pendulum" 
    annotation (Placement(transformation(extent={{100,12},{120,32}},
          rotation=0), iconTransformation(extent={{100,32},{120,52}})));
    Modelica.Blocks.Interfaces.RealOutput x3 "Horziontal position of crab" 
    annotation (Placement(transformation(extent={{100,-12},{120,8}},
          rotation=0), iconTransformation(extent={{100,-12},{120,8}})));
    Modelica.Blocks.Interfaces.RealOutput x4 "Horziontal velocity of crab" 
    annotation (Placement(transformation(extent={{100,60},{120,80}},
          rotation=0), iconTransformation(extent={{100,-56},{120,-36}})));

    Modelica.Blocks.Interfaces.RealOutput y1 "horizontal position of load" 
    annotation (Placement(transformation(extent={{100,-100},{120,-80}},
          rotation=0), iconTransformation(extent={{100,-100},{120,-80}})));

  inner Modelica.Mechanics.MultiBody.World world(
    g=g,
    animateWorld=false,
    animateGravity=false) annotation (Placement(transformation(extent={{-72,
            -40},{-52,-20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(        animation=
       false, useAxisFlange=true) 
              annotation (Placement(transformation(extent={{-32,-40},{-12,
            -20}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute(useAxisFlange=true) 
    annotation (Placement(transformation(extent={{8,-40},{28,-20}},
          rotation=0)));
  Modelica.Mechanics.Rotational.Components.Damper damper(
                                              d=d) 
    annotation (Placement(transformation(extent={{6,-20},{26,0}},  rotation=
           0)));
  Modelica.Mechanics.Translational.Components.Mass mass(
                                                    m=m_crab) 
    annotation (Placement(transformation(extent={{-38,-10},{-18,10}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.Body body(
    m=m_load,
    I_11=0,
    I_22=0,
    I_33=J) annotation (Placement(transformation(extent={{72,-40},{92,-20}},
          rotation=0)));
  Modelica.Mechanics.Translational.Sources.Force Force1 
    annotation (Placement(transformation(extent={{-78,-10},{-58,10}},
          rotation=0)));
  Modelica.Mechanics.Translational.Sensors.PositionSensor crabPosition 
    annotation (Placement(transformation(extent={{10,82},{30,102}},
                                                                  rotation=
            0)));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation translation(r={0,-l,0}) 
    annotation (Placement(transformation(extent={{40,-40},{60,-20}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape box2(
    length=0.2,
    width=0.1,
    height=0.1,
    r_shape={-0.1,0,0}) annotation (Placement(transformation(extent={{0,-66},
            {20,-46}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape wheel1(
    length=0.2,
    shapeType="cylinder",
    r_shape={-0.05,0,0},
    lengthDirection={0,0,1},
    width=0.05,
    height=0.05) annotation (Placement(transformation(extent={{-12,-98},{8,-78}},
                   rotation=0)));
  Modelica.Mechanics.Rotational.Sensors.RelAngleSensor relAngleSensor 
    annotation (Placement(transformation(extent={{6,50},{26,30}})));

  Modelica.Mechanics.Rotational.Sensors.RelSpeedSensor relSpeedSensor 
    annotation (Placement(transformation(extent={{6,20},{26,0}})));
  Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor 
    annotation (Placement(transformation(extent={{10,60},{30,80}})));

  Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition 
    annotation (Placement(transformation(extent={{36,-84},{56,-64}})));
equation
  connect(world.frame_b, prismatic.frame_a) annotation (Line(
      points={{-52,-30},{-32,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(prismatic.frame_b, revolute.frame_a) annotation (Line(
      points={{-12,-30},{8,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(damper.flange_b, revolute.axis) annotation (Line(points={{26,-10},{26,
          -20},{18,-20}},    color={0,0,0}));
  connect(revolute.support, damper.flange_a) annotation (Line(points={{12,-20},{
          6,-20},{6,-10}},     color={0,0,0}));
  connect(prismatic.axis, mass.flange_a) annotation (Line(points={{-14,-24},
          {-14,-16},{-38,-16},{-38,0}}, color={0,127,0}));
  connect(Force1.flange,   mass.flange_a) annotation (Line(points={{-58,0},
          {-38,0}}, color={0,127,0}));
  connect(force, Force1.f) annotation (Line(points={{-120,0},{-80,0}},
        color={0,0,127}));
  connect(prismatic.axis,crabPosition.flange)    annotation (Line(points={{-14,-24},
          {-14,92},{10,92}},         color={0,127,0}));
  connect(crabPosition.s,x3)  annotation (Line(points={{31,92},{72.5,92},{72.5,-2},
          {110,-2}},     color={0,0,127}));
  connect(revolute.frame_b, translation.frame_a) annotation (Line(
      points={{28,-30},{40,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(translation.frame_b, body.frame_a) annotation (Line(
      points={{60,-30},{72,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(box2.frame_a, revolute.frame_a) annotation (Line(
      points={{0,-56},{-4,-56},{-4,-30},{8,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(wheel1.frame_a, prismatic.frame_b) annotation (Line(
      points={{-12,-88},{-12,-30}},
      color={0,0,0},
      thickness=0.5));
  connect(x1, relAngleSensor.phi_rel) annotation (Line(
      points={{110,52},{16,52},{16,51}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relAngleSensor.flange_b, revolute.axis) annotation (Line(
      points={{26,40},{30,40},{30,-20},{18,-20}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(relAngleSensor.flange_a, revolute.support) annotation (Line(
      points={{6,40},{2,40},{2,-20},{12,-20}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(x1, x1) annotation (Line(
      points={{110,52},{110,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relSpeedSensor.flange_b, revolute.axis) annotation (Line(
      points={{26,10},{30,10},{30,-20},{18,-20}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(relSpeedSensor.flange_a, revolute.support) annotation (Line(
      points={{6,10},{2,10},{2,-20},{12,-20}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(relSpeedSensor.w_rel, x2) annotation (Line(
      points={{16,21},{16,22},{110,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(speedSensor.flange, prismatic.axis) annotation (Line(
      points={{10,70},{-14,70},{-14,-24}},
      color={0,127,0},
      smooth=Smooth.None));
  connect(speedSensor.v, x4) annotation (Line(
      points={{31,70},{110,70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativePosition.frame_b, translation.frame_b) annotation (Line(
      points={{56,-74},{62,-74},{62,-30},{60,-30}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  connect(relativePosition.r_rel[1], y1) annotation (Line(
      points={{46,-85.6667},{76,-85.6667},{76,-90},{110,-90}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(relativePosition.frame_a, prismatic.frame_a) annotation (Line(
      points={{36,-74},{-32,-74},{-32,-30}},
      color={95,95,95},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
            {100,100}}),       graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-146,162},{144,102}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="%name"),
        Line(points={{-90,40},{20,40}}, color={0,0,0}),
        Rectangle(
          extent={{-80,70},{-40,40}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,40},{0,-60}}, color={0,0,0}),
        Ellipse(
          extent={{-20,-40},{20,-80}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{40,106},{100,66}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="x1"),
        Text(
          extent={{40,-24},{100,-64}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="x4"),
        Text(
          extent={{42,-68},{102,-108}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="y1"),
        Text(
          extent={{40,18},{100,-22}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="x3"),
        Text(
          extent={{40,64},{100,24}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="x2")}),
                       Diagram(graphics));
end CraneMultibody;
