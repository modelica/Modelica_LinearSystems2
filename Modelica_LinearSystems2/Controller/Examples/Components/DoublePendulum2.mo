within Modelica_LinearSystems2.Controller.Examples.Components;
model DoublePendulum2 "Crane trolley system"

  extends Modelica_LinearSystems2.Controller.Templates.Internal.PlantTemplate(
                                  n=6, l=6);
  constant Real pi=Modelica.Constants.pi;
  parameter Modelica.Units.SI.Mass m_trolley=1000 "Mass of trolley";
  parameter Modelica.Units.SI.Mass m_load=4000 "Mass of load on 2nd arm";
  parameter Modelica.Units.SI.Length length=10
    "Total length of double pendulum (i.e. length of each arm = length/2)";
  parameter Modelica.Units.SI.Angle phi1_start=-40.0/180*pi
    "Initial rotation angle of 1st arm relative to trolley";
  parameter Modelica.Units.SI.Angle phi2_start=-70.0/180*pi
    "Initial rotation angle of 2nd arm relative to 1st arm";
  parameter Modelica.Units.SI.AngularVelocity w1_start=0.0
    "Initial angular velocity of 1st arm";
  parameter Modelica.Units.SI.AngularVelocity w2_start=0.0
    "Initial angular velocity of 2nd arm";

  Modelica_LinearSystems2.Utilities.Plants.DoublePendulum
    doublePendulum(
    phi1_start=phi1_start,
    phi2_start=phi2_start,
    w1_start=w1_start,
    w2_start=w2_start,
    m_trolley=m_trolley,
    m_load=m_load,
    length=length)
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  Modelica.Blocks.Routing.Multiplex6 multiplex6_1
    annotation (Placement(transformation(extent={{20,-10},{40,10}})));
equation
  connect(multiplex6_1.u1[1], doublePendulum.s) annotation (Line(
      points={{18.8,8.5},{0,8.5},{0,10},{-19,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u2[1], doublePendulum.v) annotation (Line(
      points={{18.8,5.1},{0,5.1},{0,6},{-19,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u3[1], doublePendulum.phi) annotation (Line(
      points={{18.8,1.7},{0,1.7},{0,2},{-19,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u4[1], doublePendulum.w) annotation (Line(
      points={{18.8,-1.7},{0,-1.7},{0,-2},{-19,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u5[1], doublePendulum.phi1) annotation (Line(
      points={{18.8,-5.1},{0,-5.1},{0,-6},{-19,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u6[1], doublePendulum.w1) annotation (Line(
      points={{18.8,-8.5},{0,-8.5},{0,-10},{-19,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(doublePendulum.u, u[1]) annotation (Line(
      points={{-42,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, ym) annotation (Line(
      points={{41,0},{50,0},{50,-40},{0,-40},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, y) annotation (Line(
      points={{41,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Documentation(info="<html>
<p>Model of a simple double pendulum system using <a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\">double pendulum multibody model</a>.</p>
<p>The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.doublePendulumController where it is being linearized an used as a base for linear controller design. The results are used to control the crane system as shown in <a href=\"Modelica://Modelica_LinearSystems2.Controller.Examples.DoublePendulum\">double pendulum example</a>.</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-76,10},{88,6}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Rectangle(extent={{-38,42},{6,16}}, lineColor={0,0,0}),
        Ellipse(
          extent={{-34,22},{-22,10}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255},
          lineThickness=0.5),
        Ellipse(
          extent={{-10,22},{2,10}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Line(
          points={{-12,-28},{16,-74}},
          color={0,0,0},
          smooth=Smooth.None),
        Ellipse(
          extent={{10,-68},{26,-84}},
          lineColor={0,0,0},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-19,32},{-13,26}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{34,34},{10,34}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{40,28},{16,28}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-16,28},{-12,-28}},
          color={0,0,0},
          smooth=Smooth.None),
        Ellipse(
          extent={{-14,-27},{-8,-33}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,140},{150,100}},
          lineColor={0,0,255},
          textString="%name")}));
end DoublePendulum2;
