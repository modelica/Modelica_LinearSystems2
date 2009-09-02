within Modelica_LinearSystems2.Controller;
block FirstOrder
  "First order (continuous or discrete) transfer function block (= 1 pole)"
  extends Interfaces.PartialSISO2(y(start=y_start), discretePart(
      x_start={y_start},
      y_start={y_start},
      ABCD=[-1/T,k/T; 1,0]));

  parameter Real k = 1 "Gain";
  parameter Real T = 1 "Time Constant";
  parameter Real y_start = 0.0
    "Initial y if initType=InitialState (else guess)"
    annotation(Dialog(tab="Advanced options"));

  annotation (
  defaultComponentName="firstOrder",
    Window(
      x=0.27,
      y=0.1,
      width=0.57,
      height=0.75),
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-80,80},{-80,-88}}, color={192,192,192}),
        Polygon(
          points={{-80,92},{-88,70},{-72,70},{-80,90},{-80,92}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-78},{82,-78}}, color={192,192,192}),
        Polygon(
          points={{90,-78},{68,-70},{68,-86},{90,-78}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-22,0},{88,-62}},
          lineColor={192,192,192},
          textString="PT1"),
        Text(
          extent={{-98,-106},{98,-146}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="T=%T"),
        Line(points={{-80,-78},{-70,-43.11},{-60,-17.58},{-50,1.0913},{-40,
              14.75},{-30,24.75},{-20,32.06},{-10,37.41},{0,41.33},{10,44.19},{
              20,46.29},{30,47.82},{40,48.94},{50,49.76},{60,50.36},{70,50.8},{
              80,51.12}}, color={0,0,127})}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Text(
          extent={{-48,52},{50,8}},
          lineColor={0,0,0},
          textString="k"),
        Text(
          extent={{-54,-6},{56,-56}},
          lineColor={0,0,0},
          textString="T s + 1"),
        Line(points={{-50,0},{50,0}}, color={0,0,0}),
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,127}),
        Line(points={{-100,0},{-60,0}}, color={0,0,127}),
        Line(points={{60,0},{100,0}}, color={0,0,127})}),
    Documentation(info="
<HTML>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>first order</i> system:
</p>
<pre>             k
     y = --------- * u
         T * s + 1
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
</HTML>
"));

equation
  if continuous then
    der(y) = (k*u - y)/T;
  end if;
  connect(y,discretePart.x[1]);
initial equation
  if continuous then
     if init == Types.Init.InitialState or init == Types.Init.InitialOutput then
        y = y_start;
     elseif init == Types.Init.SteadyState then
        der(y) = 0;
     end if;
  end if;
end FirstOrder;
