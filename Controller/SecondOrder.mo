within Modelica_LinearSystems2.Controller;
block SecondOrder
  "Second order (continuous or discrete) transfer function block (= 2 poles)"
  extends Interfaces.PartialSISO2(y(start=y_start), discretePart(
      x_start={y_start,yd_start},
      y_start={y_start},
      ABCD=[[0,1; -w*w,-2*D*w],[0; w*w*k]; [1,0],[0]]));
  Modelica.Blocks.Interfaces.RealOutput yd(start=yd_start)
    "First derivative of y";
  Modelica.Blocks.Interfaces.RealOutput yy[2]={y,yd};
  parameter Real k=1 "Gain";
  parameter Real w=1 "Angular frequency";
  parameter Real D=1 "Damping";
  parameter Real y_start=0.0 "Initial or guess value of output (= state)"
    annotation(Dialog(tab="Advanced options"));
  parameter Real yd_start=0.0
    "Initial or guess value of derivative of output (= state)"
    annotation(Dialog(tab="Advanced options"));

equation
  if continuous then
    der(y) = yd;
    der(yd) = w*(w*(k*u - y) - 2*D*yd);
  end if;

  connect(yy, discretePart.x);

initial equation
  if continuous then
    if init == Types.Init.SteadyState then
      der(y) = 0;
      der(yd) = 0;
    elseif init == Types.Init.InitialState then
      y = y_start;
      yd = yd_start;
    elseif init == Types.Init.InitialOutput then
      y = y_start;
      yd = 0;
    end if;
  end if;
  annotation (
    defaultComponentName="secondOrder",
    Documentation(info="
<HTML>
<p>
This blocks defines the transfer function between the input u and
the output y as <i>second order</i> system:
</p>
<pre>
                         k
     y = --------------------------------- * u
          ( s / w )^2 + 2*D*( s / w ) + 1
</pre>
<p>
The block can be continuous or discrete (with continuous parameterization).
</p>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general model class <b>TransferFunction</b>
instead and model a second order SISO system with parameters<br>
n = {k}, d = {1/w^2, 2*D/w, 1}.
</p>
<pre>
Example:
   parameter: k =  0.3,  w = 0.5,  D = 0.4
   results in:
                  0.3
      y = ------------------- * u
          4.0 s^2 + 1.6 s + 1
</pre>
</HTML>
"), Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,88},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-72,-68.53},{-64,-39.5},{-56,-2.522},{-48,32.75},
              {-40,58.8},{-32,71.51},{-24,70.49},{-16,58.45},{-8,40.06},{0,
              20.55},{8,4.459},{16,-5.271},{24,-7.629},{32,-3.428},{40,5.21},{
              48,15.56},{56,25.03},{64,31.66},{72,34.5},{80,33.61}}, color={0,0,
              127}),
        Text(
          extent={{-42,-28},{92,-68}},
          lineColor={192,192,192},
          textString="PT2"),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="w=%w"),
        Text(
          extent={{-78,92},{98,60}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}),
        Text(
          extent={{-60,60},{60,14}},
          lineColor={0,0,0},
          textString="k"),
        Text(
          extent={{-60,8},{-32,-20}},
          lineColor={0,0,0},
          textString="s"),
        Line(points={{-100,0},{-60,0}}),
        Line(points={{60,0},{100,0}}),
        Line(points={{-50,14},{50,14}}, color={0,0,0}),
        Line(points={{-54,-20},{-38,-20}}, color={0,0,0}),
        Text(
          extent={{-52,-26},{-36,-48}},
          lineColor={0,0,0},
          textString="w"),
        Line(points={{-50,2},{-56,-8},{-56,-28},{-52,-46}}, color={0,0,0}),
        Line(points={{-40,2},{-34,-10},{-34,-30},{-38,-46}}, color={0,0,0}),
        Text(
          extent={{-34,8},{-22,-10}},
          lineColor={0,0,0},
          textString="2"),
        Text(
          extent={{-34,-6},{6,-36}},
          lineColor={0,0,0},
          textString="+2D"),
        Text(
          extent={{2,8},{30,-20}},
          lineColor={0,0,0},
          textString="s"),
        Line(points={{8,-20},{24,-20}}, color={0,0,0}),
        Text(
          extent={{10,-26},{26,-48}},
          lineColor={0,0,0},
          textString="w"),
        Line(points={{12,2},{6,-8},{6,-28},{10,-46}}, color={0,0,0}),
        Line(points={{22,2},{28,-10},{28,-30},{24,-46}}, color={0,0,0}),
        Text(
          extent={{30,2},{58,-42}},
          lineColor={0,0,0},
          textString="+1")}));
end SecondOrder;
