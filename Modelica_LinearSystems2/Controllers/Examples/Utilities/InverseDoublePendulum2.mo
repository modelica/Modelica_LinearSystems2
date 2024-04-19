within Modelica_LinearSystems2.Controllers.Examples.Utilities;
model InverseDoublePendulum2 "Inverted double pendulum"
  extends Modelica_LinearSystems2.Controllers.Templates.PartialPlantMIMO(
    n=6,
    l=6);

  parameter Modelica.Units.SI.Mass m_trolley = 1 "Mass of trolley";
  parameter Modelica.Units.SI.Mass m_load = 4 "Mass of load on 2nd arm";
  parameter Modelica.Units.SI.Length length = 1
    "Total length of double pendulum (i.e. length of each arm = length/2)";
  parameter Modelica.Units.SI.Position s_start = 0.0
    "Initial position of trolley relative to world";
  parameter Modelica.Units.SI.Velocity v_start = 0.0
    "Initial velocity of trolley relative to world";
  parameter Modelica.Units.SI.Angle phi1_start = 90.0/180*pi
    "Initial rotation angle of 1st arm relative to trolley";
  parameter Modelica.Units.SI.Angle phi2_start = 0
    "Initial rotation angle of 2nd arm relative to 1st arm";
  parameter Modelica.Units.SI.AngularVelocity w1_start = 0.0
    "Initial angular velocity of 1st arm relative to trolley";
  parameter Modelica.Units.SI.AngularVelocity w2_start = 0.0
    "Initial angular velocity of 2nd arm relative to 1st arm";

  parameter Boolean cartDisturbance=false
    "True, if cart disturbance should be enabled";
  parameter Boolean bodyDisturbance=false
    "True, if body disturbance should be enabled";

  constant Real pi=Modelica.Constants.pi;

  Real dist1_s;
  Real dist2_v=0;
  Real dist3_phi1=0;
  Real dist4_w1=0;
  Real dist5_phi2=0;
  Real dist6_w2=0;

  InverseDoublePendulum inverseDoublePendulum(
    m_trolley=m_trolley,
    m_load=m_load,
    length=length,
    s_start=s_start,
    v_start=v_start,
    phi1_start=phi1_start,
    phi2_start=phi2_start,
    w1_start=w1_start,
    w2_start=w2_start,
    cartDisturbance=cartDisturbance,
    bodyDisturbance=bodyDisturbance) annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  Modelica.Blocks.Routing.Multiplex6 multiplex6
    annotation (Placement(transformation(extent={{20,-10},{40,10}})));
  Modelica.Blocks.Interfaces.RealInput dist if cartDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-60,120})));
  Modelica.Blocks.Interfaces.RealInput dist2 if bodyDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={60,120})));
  Internal.Add add[6] annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=-90,
        origin={0,-60})));
  Modelica.Blocks.Sources.RealExpression measureDisturbance[6](
    y={dist1_s,dist2_v,dist3_phi1,dist4_w1,dist5_phi2,dist6_w2})
    annotation (Placement(transformation(extent={{-40,-70},{-20,-50}})));

equation
  dist1_s=0.02*Modelica.Math.sin(20*time);
  connect(measureDisturbance.y, add.u1) annotation (Line(
      points={{-19,-60},{-8,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inverseDoublePendulum.s, multiplex6.u1[1]) annotation (Line(
      points={{-19,10},{-0.1,10},{-0.1,8.5},{18.8,8.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u2[1], inverseDoublePendulum.v) annotation (Line(
      points={{18.8,5.1},{0,5.1},{0,6},{-19,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u3[1], inverseDoublePendulum.phi) annotation (Line(
      points={{18.8,1.7},{0,1.7},{0,2},{-19,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u4[1], inverseDoublePendulum.w) annotation (Line(
      points={{18.8,-1.7},{0,-1.7},{0,-2},{-19,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u5[1], inverseDoublePendulum.phi1) annotation (Line(
      points={{18.8,-5.1},{0,-5.1},{0,-6},{-19,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u6[1], inverseDoublePendulum.w1) annotation (Line(
      points={{18.8,-8.5},{0,-8.5},{0,-10},{-19,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, ym) annotation (Line(
      points={{-1.65327e-015,-69},{-1.65327e-015,-90},{0,-90},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.y, y) annotation (Line(
      points={{41,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist2, inverseDoublePendulum.dist2) annotation (Line(
      points={{60,120},{60,60},{-24,60},{-24,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist, inverseDoublePendulum.dist) annotation (Line(
      points={{-60,120},{-60,60},{-36,60},{-36,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u[1], inverseDoublePendulum.u) annotation (Line(
      points={{-120,0},{-42,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, multiplex6.y) annotation (Line(
      points={{1.46958e-015,-52},{1.46958e-015,-40},{50,-40},{50,0},{41,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Documentation(info="<html>
<p>
Model of a&nbsp;simple inverted double pendulum system using
<a href=\"modelica://Modelica_LinearSystems2.Controllers.Examples.Utilities.InverseDoublePendulum\">InverseDoublePendulum</a>.

<!-- The following text was commented by Tobolar (DLR) since the example
cited is still in WorkInProgress (see also documentation source of
other examples of double pendulum):
The physical Model is used in
Modelica_LinearSystems2.Examples.StateSpace.inverseDoublePendulumController
where it is being linearized an used as a base for linear controller design.
-->
It is used, e.g., to control the crane system in
<a href=\"modelica://Modelica_LinearSystems2.Controllers.Examples.InverseDoublePendulum\">Controllers.Examples.InverseDoublePendulum</a>.
</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
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
end InverseDoublePendulum2;
