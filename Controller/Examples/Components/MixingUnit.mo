within Modelica_LinearSystems2.Controller.Examples.Components;
model MixingUnit
  "Mixing unit demo from Foellinger, Nichtlineare Regelungen II, p. 280"
  extends Templates.Internal.PlantTemplate_SISO;

  parameter Real a1 = 0.2674;
  parameter Real a21 = 1.815;
  parameter Real a22 = 0.4682;
  parameter Real b = 1.5476;
  parameter Real k0 = 1.05e14;
  parameter Real eps = 34.2894;
  parameter Real c0(unit="mol/l") = 0.848 "Nominal concentration";
  parameter Modelica.SIunits.Temperature T0 = 308.5 "Nominal temperature";
  parameter Real x10 = 0.42;
  parameter Real x10_inv = 0.6;
  parameter Real x20 = 0.01;
  parameter Real u0 = -0.021439;

  final parameter Real c_start(unit="mol/l") = c0*(1-x10);
  final parameter Real c_inv_start(unit="mol/l") = c0*(1-x10_inv);
  final parameter Modelica.SIunits.Temperature T_start = T0*(1+x20);
  final parameter Real c_high_start(unit="mol/l") = c0*(1-0.72);
  final parameter Real T_c_start = T0*(1+u0);

  MixingUnit1 mixingUnit(
    c0=c0,
    T0=T0,
    a1=a1,
    a21=a21,
    a22=a22,
    b=b,
    k0=k0,
    eps=eps)
    annotation (Placement(transformation(extent={{-40,-40},{40,40}})));

equation
  connect(mixingUnit.c, y) annotation (Line(
      points={{48,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(mixingUnit.T, ym[1]) annotation (Line(
      points={{0,-48},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, mixingUnit.T_c) annotation (Line(
      points={{-120,0},{-48,0}},
      color={0,0,127},
      smooth=Smooth.None));
annotation (                       Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,40},{100,-100}},
          lineColor={255,255,255},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,0}),
        Text(
          extent={{-144,148},{148,100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="%name"),
        Text(
          extent={{8,-100},{50,-138}},
          lineColor={0,0,0},
          textString="T"),
        Text(
          extent={{92,44},{134,6}},
          lineColor={0,0,0},
          textString="c"),
        Text(
          extent={{-194,72},{-108,34}},
          lineColor={0,0,0},
          textString="T_c"),
        Line(points={{0,-50},{0,-100}}, color={0,0,0}),
        Ellipse(
          extent={{-42,-38},{0,-66}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{0,-38},{42,-66}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{0,80},{0,16}}, color={0,0,0}),
        Line(points={{20,80},{20,16}}, color={0,0,0}),
        Line(points={{-86,-72},{-86,-114}}, color={0,0,0}),
        Line(points={{-66,-72},{-66,-114}}, color={0,0,0})}));
end MixingUnit;
