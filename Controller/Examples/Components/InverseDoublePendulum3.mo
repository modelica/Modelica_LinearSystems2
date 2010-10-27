within Modelica_LinearSystems2.Controller.Examples.Components;
model InverseDoublePendulum3 "Inverted double pendulum"

  extends Modelica_LinearSystems2.Controller.Templates.Internal.PlantTemplate(
                                                                     n=6, l=if secondAngle then 3 else 2);
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

  parameter Boolean secondAngle=false;

  Components.InverseDoublePendulum inverseDoublePendulum(
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
Real dist1_s;
Real dist2_v=0;
Real dist3_phi1=0;
Real dist4_w1=0;
Real dist5_phi2=0;
Real dist6_w2=0;

  Modelica.Blocks.Routing.Multiplex3 multiplex3_1 if secondAngle
    annotation (Placement(transformation(extent={{40,-40},{60,-20}})));
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
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1 if not secondAngle
    annotation (Placement(transformation(extent={{10,-66},{30,-46}})));
equation
dist1_s=0.02*Modelica.Math.sin(20*time);
  connect(multiplex6_1.u1[1], inverseDoublePendulum.s)
                                                annotation (Line(
      points={{38.8,8.5},{19.4,8.5},{19.4,10},{1,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u2[1], inverseDoublePendulum.v)
                                                annotation (Line(
      points={{38.8,5.1},{19.4,5.1},{19.4,6},{1,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u3[1], inverseDoublePendulum.phi)
                                                  annotation (Line(
      points={{38.8,1.7},{20.4,1.7},{20.4,2},{1,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u4[1], inverseDoublePendulum.w)
                                                annotation (Line(
      points={{38.8,-1.7},{20.4,-1.7},{20.4,-2},{1,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u5[1], inverseDoublePendulum.phi1)
                                                   annotation (Line(
      points={{38.8,-5.1},{19.4,-5.1},{19.4,-6},{1,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.u6[1], inverseDoublePendulum.w1)
                                                 annotation (Line(
      points={{38.8,-8.5},{20.4,-8.5},{20.4,-10},{1,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inverseDoublePendulum.u, u[1])
                                  annotation (Line(
      points={{-32,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex6_1.y, y) annotation (Line(
      points={{61,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.u1[1], inverseDoublePendulum.s) annotation (Line(
      points={{38,-23},{14,-23},{14,10},{1,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.u2[1], inverseDoublePendulum.phi) annotation (Line(
      points={{38,-30},{12,-30},{12,2},{1,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.u3[1], inverseDoublePendulum.phi1) annotation (Line(
      points={{38,-37},{10,-37},{10,-6},{1,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.y, ym) annotation (Line(
      points={{61,-30},{70,-30},{70,-80},{0,-80},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist, inverseDoublePendulum.dist) annotation (Line(
      points={{-60,80},{-60,28},{-23,28},{-23,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dist2, inverseDoublePendulum.dist2) annotation (Line(
      points={{60,80},{60,50},{-7,50},{-7,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u1[1], inverseDoublePendulum.s) annotation (Line(
      points={{8,-50},{4,-50},{4,10},{1,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2[1], inverseDoublePendulum.phi) annotation (Line(
      points={{8,-62},{6,-62},{6,2},{1,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.y, ym) annotation (Line(
      points={{31,-56},{42,-56},{42,-74},{0,-74},{0,-110}},
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
end InverseDoublePendulum3;
