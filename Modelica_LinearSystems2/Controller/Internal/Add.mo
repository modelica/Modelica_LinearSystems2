within Modelica_LinearSystems2.Controller.Internal;
model Add "Output the sum of the two real inputs"

  Modelica.Blocks.Interfaces.RealInput u1 annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={0,80})));
  Modelica.Blocks.Interfaces.RealInput u2
    annotation (Placement(transformation(extent={{-100,-20},{-60,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
equation
  y = u1+u2;
  annotation ( Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
        Ellipse(
          extent={{-20,20},{20,-20}},
          lineColor={0,0,127},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid),
        Line(
          visible=fromLeft,
          points={{-100,0},{-20,0}},
          color={0,0,127}),
        Line(points={{20,0},{100,0}}, color={0,0,127}),
        Line(points={{0,100},{0,22}}, color={0,0,127}),
        Text(
          extent={{-150,-60},{150,-100}},
          lineColor={0,0,255},
          textString="%name")}),
    Documentation(info="<html>
<p>
This blocks computes output <b>y</b> as <em>sum</em> of the
two input signals <b>u1</b> and <b>u2</b>:
</p>
<pre>
    <b>y</b> = <b>u1</b> + <b>u2</b>;
</pre>
</html>"));
end Add;
