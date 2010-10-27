within Modelica_LinearSystems2.Controller.Examples.Components;
model InverseDoublePendulum2 "inverted double pendulum"
  import Modelica_LinearSystems2.Controller;
  extends Templates.Internal.PlantTemplate(n=6, l=6);
  constant Real pi=Modelica.Constants.pi;
  parameter Modelica.SIunits.Mass m_trolley=1;
  parameter Modelica.SIunits.Mass m_load=4;
  parameter Modelica.SIunits.Length length=1;
  parameter Modelica.SIunits.Angle phi1_start=90.0/180*pi;
  parameter Modelica.SIunits.Angle phi2_start=0;
  parameter Modelica.SIunits.AngularVelocity w1_start=0.0;
  parameter Modelica.SIunits.AngularVelocity w2_start=0.0;
  parameter Modelica.SIunits.Position s_start=0.0;

  parameter Boolean cartDisturbance=false;
  parameter Boolean bodyDisturbance=false;

  Controller.Examples.Components.InverseDoublePendulum inverseDoublePendulum(
    s_start=s_start,
    phi1_start=phi1_start,
    phi2_start=phi2_start,
    w1_start=w1_start,
    w2_start=w2_start,
    m_trolley=m_trolley,
    m_load=m_load,
    length=length,
    cartDisturbance=cartDisturbance,
    bodyDisturbance=bodyDisturbance)
    annotation (Placement(transformation(extent={{-30,-10},{0,10}})));
  Modelica.Blocks.Routing.Multiplex6 multiplex6_1
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
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
        origin={40,-50})));
  Modelica.Blocks.Sources.RealExpression measureDisturbance[6](y={dist1_s,
        dist2_v,dist3_phi1,dist4_w1,dist5_phi2,dist6_w2})
    annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
Real dist1_s;
Real dist2_v=0;
Real dist3_phi1=0;
Real dist4_w1=0;
Real dist5_phi2=0;
Real dist6_w2=0;

equation
dist1_s=0.02*Modelica.Math.sin(20*time);
  connect(measureDisturbance.y, add.u1) annotation (Line(
      points={{-19,-50},{32,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inverseDoublePendulum.s, multiplex6_1.u1[1]) annotation (Line(
      points={{1,10},{19.9,10},{19.9,8.5},{38.8,8.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u2[1], inverseDoublePendulum.v) annotation (Line(
      points={{38.8,5.1},{19.4,5.1},{19.4,6},{1,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u3[1], inverseDoublePendulum.phi) annotation (Line(
      points={{38.8,1.7},{20.4,1.7},{20.4,2},{1,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u4[1], inverseDoublePendulum.w) annotation (Line(
      points={{38.8,-1.7},{20.4,-1.7},{20.4,-2},{1,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u5[1], inverseDoublePendulum.phi1) annotation (Line(
      points={{38.8,-5.1},{20.4,-5.1},{20.4,-6},{1,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u6[1], inverseDoublePendulum.w1) annotation (Line(
      points={{38.8,-8.5},{20.4,-8.5},{20.4,-10},{1,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, ym) annotation (Line(
      points={{40,-59},{40,-80},{0,-80},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, y) annotation (Line(
      points={{61,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist2, inverseDoublePendulum.dist2) annotation (Line(
      points={{60,120},{60,60},{-7,60},{-7,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist, inverseDoublePendulum.dist) annotation (Line(
      points={{-60,120},{-60,60},{-23,60},{-23,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u[1], inverseDoublePendulum.u) annotation (Line(
      points={{-120,0},{-32,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, multiplex6_1.y) annotation (Line(
      points={{40,-42},{40,-28},{76,-28},{76,0},{61,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=20),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
 
 
Model of a simple inverted double pendulum system using Modelica_Controller.Examples.Components.InverseDoublePendulum.<br>
The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.inverseDoublePendulumController where it is being
linearized an used as a base for linear controller design. The results are used to control the crane system
in Modelica_Controller.Examples.InverseDoublePendulum.mo
</html>"),
    uses(Modelica(version="3.0")),
    experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{150,
            100}}), graphics={
        Rectangle(
          extent={{-150,100},{150,-100}},
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
