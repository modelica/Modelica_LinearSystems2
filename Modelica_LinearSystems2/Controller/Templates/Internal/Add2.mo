within Modelica_LinearSystems2.Controller.Templates.Internal;
block Add2
  "Output the sum of the two inputs (inputs are on the left and below)"

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
This blocks computes output <strong>y</strong> as <em>difference</em> of the
commanded input <strong>u1</strong> and the feedback
input <strong>u2</strong>:
</p>
<pre>
    <strong>y</strong> = <strong>u1</strong> - <strong>u2</strong>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   n = 2

  results in the following equations:

     y = u1 - u2
</pre>
</html>"), Icon(coordinateSystem(
    preserveAspectRatio=true,
    extent={{-100,-100},{100,100}},
    grid={2,2}), graphics={
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
          lineColor={0,0,255})}),
    Diagram(coordinateSystem(
    preserveAspectRatio=true,
    extent={{-100,-100},{100,100}},
    grid={2,2}), graphics={
        Ellipse(
          extent={{-20,20},{20,-20}},
          pattern=LinePattern.Solid,
          lineThickness=0.25,
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,255}),
        Line(points={{-60,0},{-20,0}}, color={0,0,255}),
        Line(points={{20,0},{80,0}}, color={0,0,255}),
        Line(points={{0,-20},{0,-60}}, color={0,0,255})}));
end Add2;
