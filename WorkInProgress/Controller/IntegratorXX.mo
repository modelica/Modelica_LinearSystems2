within Modelica_LinearSystems2.WorkInProgress.Controller;
block IntegratorXX
  "Output the integral of the input signal (continuous or discrete block)"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Controller.Types;
  extends Modelica_LinearSystems2.Controller.Interfaces.PartialSISO2(y(start=
          y_start), discretePart(
      withDelay=withDelay,
      x_start={y_start},
      y_start={y_start},
      ABCD=[0,k; 1,0]));

  parameter Real k=1 "Integrator gain";
  parameter Boolean withDelay=false
    "= true, if the output is delayed by one sample period (only if discrete)";
  parameter Boolean limitsAtInit=true
    "= false, if limits are ignored during initializiation (i.e., y=u)";

  parameter Real y_start=0 "Initial or guess value of output (=state)"
                                                               annotation(Dialog(tab="Advanced options"));
  Modelica.Blocks.Interfaces.RealInput lowerLimit "lower limit"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,-80})));
  Modelica.Blocks.Interfaces.RealInput upperLimit "upper limit"
    annotation (extent=[-140, -20; -100, 20], Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,80})));
protected
  Boolean reset;
  Real limit;
  Modelica.Blocks.Interfaces.RealOutput y2;
equation
  if continuous then
    der(y) = if y < lowerLimit or y > upperLimit then 0 else k*u;
    der(y2) = 0;
  else
    der(y2) = if y2 < lowerLimit or y2 > upperLimit then 0 else k*u;
    der(y) = 0;
  end if;
  reset = y < lowerLimit or y > upperLimit;
  limit = if y < lowerLimit then lowerLimit else if y > upperLimit then
    upperLimit else y;
  when reset then
    reinit(y, limit);
    reinit(y2, limit);
  end when;
  if not continuous then
    y=y2;
  end if;

  connect(y2, discretePart.y[1]);
//  connect(y2, y);
initial equation
  if continuous then
    if init == Types.Init.InitialState or init == Types.Init.InitialOutput then
      y = if limitsAtInit then min(upperLimit, max(y_start, lowerLimit)) else
        y_start;
    elseif init == Types.Init.SteadyState then
      der(y) = 0;
    end if;
  end if;
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
        preserveAspectRatio=true,
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
end IntegratorXX;
