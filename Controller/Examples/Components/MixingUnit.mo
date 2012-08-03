within Modelica_LinearSystems2.Controller.Examples.Components;
model MixingUnit
  "Stirred tank reactor as controller plant with standard interfaces"
  extends Templates.Internal.PlantTemplate_SISO;
  import SI = Modelica.SIunits;

  parameter Real c0(unit="mol/l") = 0.848 "Nominal concentration"; // Should be Modelica.SIunits.Concentration in [mol/m^3]
  parameter SI.Temperature T0 = 308.5 "Nominal temperature";
  parameter Real a1 = 0.2674 "Matrix A: Coefficient of element [1,1]";
  parameter Real a21 = 1.815 "Matrix A: Coefficient of element [2,1]";
  parameter Real a22 = 0.4682 "Matrix A: Coefficient of element [2,2]";
  parameter Real b = 1.5476 "Matrix B: Coefficient of element [2,1]";
  parameter Real k0(min=0) = 1.05e14 "Constant k0 for gamma";
  parameter Real eps(min=0) = 34.2894 "Constant eps for gamma";
  parameter Real x10 = 0.42
    "Initial value of state x1 (related concentration of substance A in tank)";
  parameter Real x10_inv = 0.6 "Initial value of state x1 of inverted model";
  parameter Real x20 = 0.01
    "Initial value of state x2 (related temperature in tank)";
  parameter Real u0 = -0.021439
    "Initial related temperature of cooling medium [-]";

  final parameter Real c_start(unit="mol/l") = c0*(1-x10)
    "Initial concentration of substance A in tank";
  final parameter Real c_inv_start(unit="mol/l") = c0*(1-x10_inv)
    "Initial concentration of substance A in tank";
  final parameter SI.Temperature T_start = T0*(1+x20)
    "Initial temperature in tank";
  //final parameter Real c_high_start(unit="mol/l") = c0*(1-0.72);
  final parameter SI.Temperature T_c_start = T0*(1+u0)
    "Initial temperature of cooling medium";

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
      points={{44,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(mixingUnit.T, ym[1]) annotation (Line(
      points={{0,-44},{0,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, mixingUnit.T_c) annotation (Line(
      points={{-120,0},{-48,0}},
      color={0,0,127},
      smooth=Smooth.None));

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false,
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
          lineColor={95,95,95},
          textString="T"),
        Text(
          extent={{100,50},{142,12}},
          lineColor={95,95,95},
          textString="c"),
        Text(
          extent={{-186,68},{-100,30}},
          lineColor={95,95,95},
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
        Line(points={{-60,80},{-60,16}},
                                     color={0,0,0}),
        Line(points={{-40,80},{-40,16}},
                                       color={0,0,0}),
        Line(points={{60,-78},{60,-120}},   color={0,0,0}),
        Line(points={{80,-78},{80,-120}},   color={0,0,0})}), Documentation(
        info="<html>
<p>Model of idealized stirred tank reactor, see <a href=\"modelica://Modelica_LinearSystems2.Controller.Examples.Components.MixingUnit1\">MixingUnit1</a>
for more details.
It is intended for replacement of a plant instance in general model of controller.
</p>
</html>"));
end MixingUnit;
