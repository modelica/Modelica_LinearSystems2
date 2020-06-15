within Modelica_LinearSystems2.Controller.Examples.Components;
model MixingUnit1 "Stirred tank reactor"
  import Modelica.Units.SI;
// stirred tank reactor // stirrer vessel reactor

  parameter Real c0(unit="mol/l") = 0.848 "Nominal concentration"; // Should be Modelica.Units.SI.Concentration in [mol/m^3]
  parameter SI.Temperature T0 = 308.5 "Nominal temperature";
  parameter Real a1 = 0.2674 "Matrix A: Coefficient of element [1,1]";
  parameter Real a21 = 1.815 "Matrix A: Coefficient of element [2,1]";
  parameter Real a22 = 0.4682 "Matrix A: Coefficient of element [2,2]";
  parameter Real b = 1.5476 "Matrix B: Coefficient of element [2,1]";
  parameter Real k0(min=0) = 1.05e14 "Constant k0 for gamma";
  parameter Real eps(min=0) = 34.2894 "Constant eps for gamma";
  Real gamma "Reaction speed";
protected
  parameter SI.Time tau0 = 60 "Time scaling";
  parameter Real wk0 = k0/c0;
  parameter Real weps = eps*T0;
  parameter Real wa11 = a1/tau0;
  parameter Real wa12 = c0/tau0;
  parameter Real wa13 = c0*a1/tau0;
  parameter Real wa21 = a21/tau0;
  parameter Real wa22 = a22*T0/tau0;
  parameter Real wa23 = T0*(a21 - b)/tau0;
  parameter Real wb = b/tau0;

public
  Modelica.Blocks.Interfaces.RealInput T_c(unit="K") "Cooling temperature"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput c(unit="mol/l")
    "Concentration in mixing unit"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput T(unit="K")
    "Temperature in mixing unit"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=-90,
        origin={0,-110})));
equation
  gamma = c*wk0*exp( -weps/(T+1));
  der(c) = -wa11*c - wa12*gamma + wa13;
  der(T) = -wa21*T + wa22*gamma + wa23 + wb*T_c;
  annotation (                       Icon(coordinateSystem(preserveAspectRatio=true,
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
        Line(points={{40,80},{40,16}}, color={0,0,0}),
        Line(points={{60,80},{60,16}}, color={0,0,0}),
        Line(points={{-80,-78},{-80,-120}}, color={0,0,0}),
        Line(points={{-60,-78},{-60,-120}}, color={0,0,0}),
        Text(
          extent={{-186,68},{-100,30}},
          lineColor={95,95,95},
          textString="T_c"),
        Text(
          extent={{100,50},{142,12}},
          lineColor={95,95,95},
          textString="c"),
        Text(
          extent={{8,-100},{50,-138}},
          lineColor={95,95,95},
          textString="T")}),
    Diagram(graphics={
        Rectangle(
          extent={{-70,16},{70,-60}},
          lineColor={255,255,255},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,20},{60,-60}},
          lineColor={255,255,255},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,40},{-60,-60},{60,-60},{60,40}},         color={0,0,0}),
        Line(points={{0,-36},{0,-66}},  color={0,0,0}),
        Ellipse(
          extent={{0,-28},{42,-44}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{20,68},{20,34}},
                                     color={0,0,0}),
        Line(points={{40,68},{40,34}}, color={0,0,0}),
        Line(points={{-50,-54},{-50,-80}},  color={0,0,0}),
        Ellipse(
          extent={{-42,-28},{0,-44}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-30,-54},{-30,-80}},  color={0,0,0}),
        Line(
          points={{-64,0},{-100,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Ellipse(
          extent={{-66,2},{-62,-2}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{100,0},{46,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Ellipse(
          extent={{44,2},{48,-2}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{20,-52},{0,-90},{0,-100}},
          color={0,0,127},
          smooth=Smooth.None),
        Ellipse(
          extent={{18,-50},{22,-54}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{30,58},{30,26}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{26,34},{30,26},{34,34}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-40,-66},{-40,-94}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-44,-86},{-40,-94},{-36,-86}},
          color={0,0,0},
          smooth=Smooth.None),
        Text(
          extent={{20,94},{40,74}},
          lineColor={0,0,0},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid,
          textString="A"),
        Text(
          extent={{36,68},{86,56}},
          lineColor={0,0,0},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid,
          textString="c0, T0"),
        Rectangle(extent={{-70,40},{70,-60}}, lineColor={0,0,0})}),
    Documentation(info="<html>
<p>
Model of idealized stirred tank reactor from [1], p. 280. In the stirred tank there is
a continuous intake of substance A which reacts with a catalyzer and brakes down into
other substances.
Since the content of the tank reactor is stired continuously, it can be supposed
a homogenous distribution of all substances within the tank.
Therefore, the tank reactor system can be described by means of concentrated
parameters.
</p>
<p>
The chemical reaction in the tank reactor is exothermic, i.e. it produces themal
energy which affects the stability of the reaction.
To stabilize it, the tank reactor is cooled. The temperature&nbsp;<code>T_c</code>
of cooling medium is the variable to be actuated by controller.
</p>
<p>
The speed &gamma; of the reaction depends on the concentration&nbsp;<code>c</code> and
temperature&nbsp;<code>T</code> as follows:
</p>
<blockquote><pre>
gamma = (1-x1) * k0 * e^(-eps/(1+x2))
</pre></blockquote>
<p>
with related concentration&nbsp;<code>x1</code> and related temperature&nbsp;<code>x2</code>:
</p>
<blockquote><pre>
      c0 - c            T - T0
x1 = --------;    x2 = --------.
        c0                T0
</pre></blockquote>
<p>
Both variables <code>x1</code> and &nbsp;<code>x2</code> are states of the model.
Additionally, there is used related temperature&nbsp;<code>u</code> of cooling medium:
</p>
<blockquote><pre>
     T_c - T0
u = --------.
       T0
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] F&ouml;llinger O. (1998):</dt>
<dd> <b>Nichtlineare Regelungen I</b>.
     8th Edition, Oldenbourg Verlag M&uuml;nchen.<br>&nbsp;</dd>
</dl>
</html>"));
end MixingUnit1;
