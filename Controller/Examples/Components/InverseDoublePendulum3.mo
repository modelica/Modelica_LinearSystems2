within Modelica_LinearSystems2.Controller.Examples.Components;
model InverseDoublePendulum3 "Inverted double pendulum"
  extends Modelica_LinearSystems2.Controller.Templates.Internal.PlantTemplate(
    n=6,
    l=if secondAngle then 3 else 2);

  parameter Modelica.SIunits.Mass m_trolley = 1 "Mass of trolley";
  parameter Modelica.SIunits.Mass m_load = 4 "Mass of load on 2nd arm";
  parameter Modelica.SIunits.Length length = 1
    "Total length of double pendulum (i.e. length of each arm = length/2)";
  parameter Modelica.SIunits.Position s_start = 0.0
    "Initial position of trolley relative to world";
  parameter Modelica.SIunits.Velocity v_start = 0.0
    "Initial velocity of trolley relative to world";
  parameter Modelica.SIunits.Angle phi1_start=90.0/180*pi
    "Initial rotation angle of 1st arm relative to trolley";
  parameter Modelica.SIunits.Angle phi2_start = 0
    "Initial rotation angle of 2nd arm relative to 1st arm";
  parameter Modelica.SIunits.AngularVelocity w1_start = 0.0
    "Initial angular velocity of 1st arm relative to trolley";
  parameter Modelica.SIunits.AngularVelocity w2_start = 0.0
    "Initial angular velocity of 2nd arm relative to 1st arm";

  parameter Boolean cartDisturbance=false
    "True, if cart disturbance should be enabled";
  parameter Boolean bodyDisturbance=false
    "True, if body disturbance should be enabled";

  parameter Boolean secondAngle=false
    "Include 2nd angle of double pendulum in output ym";

  constant Real pi=Modelica.Constants.pi;

  Real dist1_s;
  Real dist2_v=0;
  Real dist3_phi1=0;
  Real dist4_w1=0;
  Real dist5_phi2=0;
  Real dist6_w2=0;

  Components.InverseDoublePendulum inverseDoublePendulum(
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
    bodyDisturbance=bodyDisturbance)
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  Modelica.Blocks.Routing.Multiplex6 multiplex6
    annotation (Placement(transformation(extent={{20,-10},{40,10}})));

  Modelica.Blocks.Routing.Multiplex3 multiplex3 if   secondAngle
    annotation (Placement(transformation(extent={{0,-40},{20,-20}})));
  Modelica.Blocks.Interfaces.RealInput dist if cartDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-60,80}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-48,56})));
  Modelica.Blocks.Interfaces.RealInput dist2 if bodyDisturbance
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={60,80}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={54,56})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2 if   not secondAngle
    annotation (Placement(transformation(extent={{0,-80},{20,-60}})));
equation
dist1_s=0.02*Modelica.Math.sin(20*time);
  connect(multiplex6.u1[1], inverseDoublePendulum.s)
                                                annotation (Line(
      points={{18.8,8.5},{0,8.5},{0,10},{-39,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u2[1], inverseDoublePendulum.v)
                                                annotation (Line(
      points={{18.8,5.1},{0,5.1},{0,6},{-39,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u3[1], inverseDoublePendulum.phi)
                                                  annotation (Line(
      points={{18.8,1.7},{0,1.7},{0,2},{-39,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u4[1], inverseDoublePendulum.w)
                                                annotation (Line(
      points={{18.8,-1.7},{0,-1.7},{0,-2},{-39,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u5[1], inverseDoublePendulum.phi1)
                                                   annotation (Line(
      points={{18.8,-5.1},{0,-5.1},{0,-6},{-39,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.u6[1], inverseDoublePendulum.w1)
                                                 annotation (Line(
      points={{18.8,-8.5},{0,-8.5},{0,-10},{-39,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inverseDoublePendulum.u, u[1])
                                  annotation (Line(
      points={{-62,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6.y, y)   annotation (Line(
      points={{41,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3.u1[1], inverseDoublePendulum.s)   annotation (Line(
      points={{-2,-23},{-18,-23},{-18,10},{-39,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3.u2[1], inverseDoublePendulum.phi)   annotation (Line(
      points={{-2,-30},{-24,-30},{-24,2},{-39,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3.u3[1], inverseDoublePendulum.phi1)   annotation (Line(
      points={{-2,-37},{-30,-37},{-30,-6},{-39,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3.y, ym)   annotation (Line(
      points={{21,-30},{30,-30},{30,-90},{0,-90},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist, inverseDoublePendulum.dist) annotation (Line(
      points={{-60,80},{-60,40},{-56,40},{-56,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist2, inverseDoublePendulum.dist2) annotation (Line(
      points={{60,80},{60,40},{-44,40},{-44,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2.u1[1], inverseDoublePendulum.s)   annotation (Line(
      points={{-2,-64},{-18,-64},{-18,10},{-39,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2.u2[1], inverseDoublePendulum.phi)   annotation (Line(
      points={{-2,-76},{-24,-76},{-24,2},{-39,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2.y, ym)   annotation (Line(
      points={{21,-70},{30,-70},{30,-90},{0,-90},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(
          extent={{-10,-54},{80,-84}},
          lineColor={95,95,95},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.Dash),
        Rectangle(
          extent={{-10,-14},{80,-44}},
          lineColor={95,95,95},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.Dash),
        Text(
          extent={{32,-16},{78,-22}},
          lineColor={0,128,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="if secondAngle=true"),
        Text(
          extent={{32,-56},{80,-62}},
          lineColor={255,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="if secondAngle=false")}),
    Documentation(info="<html>
<p>Model of a simple inverted double pendulum system using Modelica_Controller.Examples.Components.InverseDoublePendulum. The physical Model is used in Modelica_LinearSystems2.Examples.StateSpace.inverseDoublePendulumController where it is being linearized an used as a base for linear controller design. The results are used to control the crane system in Modelica_Controller.Examples.InverseDoublePendulum.mo </p>
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
end InverseDoublePendulum3;
