within Modelica_LinearSystems2.Controller.Examples.Components;
model DoublePendulum2 "crane trolley system"

  extends Templates.PlantTemplate(n=6, l=6);
  constant Real pi=Modelica.Constants.pi;
  parameter Modelica.SIunits.Mass m_trolley=1000;
  parameter Modelica.SIunits.Mass m_load=4000;
  parameter Modelica.SIunits.Length length=10;
  parameter Modelica.SIunits.Angle phi1_start=-40.0/180*pi;
  parameter Modelica.SIunits.Angle phi2_start=-70.0/180*pi;
  parameter Modelica.SIunits.AngularVelocity w1_start=0.0;
  parameter Modelica.SIunits.AngularVelocity w2_start=0.0;

  DoublePendulum doublePendulum(
    phi1_start=phi1_start,
    phi2_start=phi2_start,
    w1_start=w1_start,
    w2_start=w2_start,
    m_trolley=m_trolley,
    m_load=m_load,
    length=length) 
    annotation (Placement(transformation(extent={{-30,-10},{0,10}})));
  Modelica.Blocks.Routing.Multiplex6 multiplex6_1 
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
equation
  connect(multiplex6_1.u1[1], doublePendulum.s) annotation (Line(
      points={{38.8,8.5},{19.4,8.5},{19.4,10},{1,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u2[1], doublePendulum.v) annotation (Line(
      points={{38.8,5.1},{19.4,5.1},{19.4,6},{1,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u3[1], doublePendulum.phi) annotation (Line(
      points={{38.8,1.7},{20.4,1.7},{20.4,2},{1,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u4[1], doublePendulum.w) annotation (Line(
      points={{38.8,-1.7},{20.4,-1.7},{20.4,-2},{1,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u5[1], doublePendulum.phi1) annotation (Line(
      points={{38.8,-5.1},{19.4,-5.1},{19.4,-6},{1,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u6[1], doublePendulum.w1) annotation (Line(
      points={{38.8,-8.5},{20.4,-8.5},{20.4,-10},{1,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(doublePendulum.u, u[1]) annotation (Line(
      points={{-32,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, ym) annotation (Line(
      points={{61,0},{66,0},{66,-80},{0,-80},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, y) annotation (Line(
      points={{61,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=20),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
 
Model of a simple double pendulum system using Modelica_LinearSystems2.Controller.Examples.Components.DoublePendulum.<br>
The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.doublePendulumController where it is being
linearized an used as a base for linear controller design. The results are used to control the crane system
in Modelica_LinearSystems2.Controller.Examples.DoublePendulum.mo
 
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
          fillPattern=FillPattern.Solid)}));
end DoublePendulum2;
