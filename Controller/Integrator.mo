within Modelica_LinearSystems2.Controller;
block Integrator
  "Output the integral of the input signal (continuous or discrete block)"
  extends Interfaces.PartialSISO2(y(start=y_start), discretePart(
      withDelay=withDelay,
      x_start={y_start},
      y_start={y_start},
      ABCD=[0,k; 1,0]));
  parameter Real k=1 "Integrator gain";
  parameter Boolean withDelay=false
    "True, if the output is delayed by one sample period (only if discrete)";

  parameter Real y_start=0 "Initial or guess value of output (=state)"
                                                               annotation(Dialog(tab="Advanced options"));
equation
  if continuous then
    der(y) = k*u;

  end if;
  connect(y, discretePart.y[1]);
initial equation
  if continuous then
    if init == Types.Init.InitialState or init == Types.Init.InitialOutput then
      y = y_start;
    elseif init == Types.Init.SteadyState then
      der(y) = 0;
    end if;
  end if;
  annotation (
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>integrator</i>:
</p>
<pre>
          k
     y = --- * u
          s
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
<p>
It is not possible to initalize a continuous integrator in steady state.
For this reason, option \"initType = SteadyState\" is ignored for
a continuous integrator block and
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
          extent={{0,-10},{60,-70}},
          lineColor={192,192,192},
          textString="I"),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="k=%k"),
        Line(points={{-80,-80},{80,80}}, color={0,0,127}),
        Text(
          extent={{-94,78},{88,46}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,127}),
        Line(points={{-100,0},{-60,0}}, color={0,0,127}),
        Line(points={{60,0},{100,0}}, color={0,0,127}),
        Text(
          extent={{-36,60},{32,2}},
          lineColor={0,0,0},
          textString="k"),
        Text(
          extent={{-32,0},{36,-58}},
          lineColor={0,0,0},
          textString="s"),
        Line(points={{-46,0},{46,0}}, color={0,0,0})}));
end Integrator;
