within Modelica_LinearSystems2.Controllers.Internal;
model Add2 "Output the sum of the two real inputs (graphics can be changed)"

  parameter Boolean fromLeft = true
    "True, if second input is left (else below)"
    annotation(choices(checkBox=true));
  parameter Integer n(min=1)=1 "Number of inputs = number of outputs";

  Modelica.Blocks.Interfaces.RealInput u1[n]
    annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={0,80})));
  Modelica.Blocks.Interfaces.RealInput u2[n] if fromLeft
    annotation (Placement(transformation(extent={{-100,-20},{-60,20}})));
  Modelica.Blocks.Interfaces.RealOutput y[n]
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  Modelica.Blocks.Interfaces.RealInput u2b[n] if not fromLeft
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-80})));
  Add add[n] annotation (Placement(transformation(extent={{10,-10},{30,10}})));
equation
//y = if fromLeft then u1+u2 else u1+u2b;

  connect(add.u2, u2b) annotation (Line(
      points={{12,0},{0,0},{0,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, u2) annotation (Line(
      points={{12,0},{-80,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, u1) annotation (Line(
      points={{20,8},{20,40},{0,40},{0,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, y) annotation (Line(
      points={{29,0},{56,0},{56,0},{90,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation ( Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
        Line(
          visible=fromLeft,
          points={{-100,0},{-20,0}},
          color={0,0,127}),
        Line(points={{20,0},{100,0}}, color={0,0,127}),
        Line(points={{0,100},{0,22}}, color={0,0,127}),
        Line(
          visible=not fromLeft,
          points={{0,-18},{0,-60}},
          color={0,0,127}),
        Text(
          extent={{-150,-20},{150,-60}},
          lineColor={0,0,255},
          textString="%name"),
        Ellipse(
          extent={{-20,20},{20,-20}},
          lineColor={0,0,127},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<html>
<p>
This blocks computes output <b>y</b> as <i>sum</i> of the
two input signals <b>u1</b> and <b>u2</b>:
</p>
<pre>
    <b>y</b> = <b>u1</b> + <b>u2</b>;
</pre>
<p>
The second input can be obtained either from left side of block or from bottom.
The corresponding connectors are <b>u2</b> or <b>u2b</b>, respectively.
</p>
</html>"));
end Add2;
