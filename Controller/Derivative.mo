within Modelica_LinearSystems2.Controller;
block Derivative "Approximate derivative (continuous or discrete block)"

  extends Interfaces.PartialSISO2(discretePart(
      y(start={y_start}),
      x_start={x_start},
      y_start={y_start},
      ABCD=[-1/T,1/T; -k/T,k/T]));

  parameter Real k=1 "Gain";
  parameter Modelica.SIunits.Time T(min=Modelica.Constants.eps) = 0.01
    "Time Constant (T>0 required; T=0 is ideal derivative block)";
  parameter Real x_start=0 "Initial or guess value of state"  annotation(Dialog(tab="Advanced options"));
  parameter Real y_start=0 "Initial or guess value of output"  annotation(Dialog(tab="Advanced options"));
  Modelica.Blocks.Interfaces.RealOutput x(start=x_start)
    "State of approximative derivative";
   annotation (
    Window(
      x=0.29,
      y=0.05,
      width=0.53,
      height=0.54),
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>approximative derivative (DT1)</i>:
</p>
<pre>
             k * s
     y = ------------ * u
            T * s + 1
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
<p>
If k=0, the state space realization of the block is
specially constructed, in order that
the D-part of PID controllers can be set to zero without
introducing numerical problems.
</p>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general model class <b>TransferFunction</b>
instead and model a DT1 system with parameters<br>
n = {k,0}, d = {T,1}.
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
        Line(points={{-80,-80},{-80,60},{-70,17.95},{-60,-11.46},{-50,-32.05},{
              -40,-46.45},{-30,-56.53},{-20,-63.58},{-10,-68.51},{0,-71.96},{10,
              -74.37},{20,-76.06},{30,-77.25},{40,-78.07},{50,-78.65},{60,-79.06}},
            color={0,0,127}),
        Text(
          extent={{0,0},{60,60}},
          lineColor={192,192,192},
          textString="DT1")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,127}),
        Line(points={{-100,0},{-60,0}}, color={0,0,127}),
        Line(points={{60,0},{100,0}}, color={0,0,127}),
        Text(
          extent={{-54,52},{50,10}},
          lineColor={0,0,0},
          textString="k s"),
        Text(
          extent={{-54,-6},{52,-52}},
          lineColor={0,0,0},
          textString="T s + 1"),
        Line(points={{-50,0},{50,0}}, color={0,0,0})}));
protected
  parameter Boolean zeroGain = abs(k) < Modelica.Constants.eps
    "= true, if k is considered to be zero";
equation
  if continuous then
     der(x) = if zeroGain then 0 else (u - x)/T;
          y = if zeroGain then 0 else (k/T)*(u - x);
else

connect(x,discretePart.x[1]);
connect(y,discretePart.y[1]);
  end if;

initial equation
  if continuous then
     if init == Types.Init.InitialState or zeroGain then
        x = x_start;
     elseif init == Types.Init.SteadyState then
        der(x) = 0;
  elseif init == Types.Init.InitialOutput then
    if zeroGain then
       x = u;
    else
       y = y_start;
    end if;
     end if;
  end if;
end Derivative;
