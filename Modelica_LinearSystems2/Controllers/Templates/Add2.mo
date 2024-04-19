within Modelica_LinearSystems2.Controllers.Templates;
block Add2 "Output the sum of the two inputs"

  input Modelica.Blocks.Interfaces.RealInput u1
    annotation (Placement(transformation(
      extent={{-100,-20},{-60,20}}, rotation=0)));
  input Modelica.Blocks.Interfaces.RealInput u2
    annotation (Placement(transformation(
    origin={0,-80},
    extent={{-20,-20},{20,20}},
    rotation=90)));
  output Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(
      extent={{80,-10},{100,10}}, rotation=0)));

equation
  y = u1 + u2;
  annotation (
    Documentation(info="<html>
<p>
This block computes output&nbsp;<strong>y</strong> as <em>sum</em> of the
inputs <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
<strong>y</strong> = <strong>u1</strong> + <strong>u2</strong>;
</pre></blockquote>
</html>"),
    Icon(
      coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Ellipse(
          extent={{-20,20},{20,-20}},
          lineColor={0,0,127},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,0},{-20,0}}, color={0,0,127}),
        Line(points={{20,0},{80,0}}, color={0,0,127}),
        Line(points={{0,-20},{0,-60}}, color={0,0,127}),
        Text(
          extent={{-150,94},{150,44}},
          textString="%name",
          textColor={0,0,255})}),
    Diagram(
      coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Ellipse(
          extent={{-20,20},{20,-20}},
          pattern=LinePattern.Solid,
          lineThickness=0.25,
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,127}),
        Line(points={{-60,0},{-20,0}}, color={0,0,127}),
        Line(points={{20,0},{80,0}}, color={0,0,127}),
        Line(points={{0,-20},{0,-60}}, color={0,0,127})}));
end Add2;
