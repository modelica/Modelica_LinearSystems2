within Modelica_LinearSystems2.Controller;
block PI "Proportional-Integral controller (continuous or discrete block)"
  extends Interfaces.PartialSISO2(discretePart(
      x_start={x_start},
      y_start={y_start},
      ABCD=[0,1/T; k,k]));

  parameter Real k=1 "Gain";
  parameter Modelica.SIunits.Time T(min=Modelica.Constants.eps)=1
    "Time Constant (T>0 required)";
  parameter Real x_start=0 "Initial or guess value of state"  annotation(Dialog(tab="Advanced options"));
parameter Real y_start=0 "Initial value of output"  annotation(Dialog(tab="Advanced options"));
  Modelica.Blocks.Interfaces.RealOutput x(start=x_start) "State of block";
   annotation (
    Window(
      x=0.29,
      y=0.05,
      width=0.53,
      height=0.54),
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>integrator</i>:
</p>
<pre>                     1
      y = k * (1 + ------ ) * u
                    T*s
               T*s + 1
        = k * --------- * u
                 T*s
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
<p>
It is not possible to initalize a continuous integrator in steady state.
For this reason, option \"initType = SteadyState\" is ignored for
a continuous PI block and
interpreted as \"initType = InitialState\".
</p>
</html>
"), Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="k=%k"),
        Line(points={{-80,-80},{-80,-20},{60,80}}, color={0,0,127}),
        Text(
          extent={{0,6},{60,-56}},
          lineColor={192,192,192},
          textString="PI"),
        Text(
          extent={{-148,-188},{152,-148}},
          lineColor={0,0,0},
          textString="T=%T")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,127}),
        Line(points={{-100,0},{-60,0}}, color={0,0,127}),
        Line(points={{60,0},{100,0}}, color={0,0,127}),
        Text(
          extent={{-70,24},{-26,-18}},
          lineColor={0,0,0},
          textString="k"),
        Text(
          extent={{-34,48},{58,0}},
          lineColor={0,0,0},
          textString="T s + 1"),
        Text(
          extent={{-32,-8},{50,-40}},
          lineColor={0,0,0},
          textString="T s"),
        Line(points={{-26,0},{52,0}}, color={0,0,0})}));
equation
  connect(y, discretePart.y[1]);
  connect(x, discretePart.x[1]);

  if continuous then
    der(x) = u/T;
    y = k*(x + u);
  end if;

initial equation
  if continuous then
    if init == Types.Init.SteadyState then
    der(x) = 0;
  elseif init == Types.Init.InitialState then
    x = x_start;
  elseif init == Types.Init.InitialOutput then
    y = y_start;
  end if;
  end if;
end PI;
