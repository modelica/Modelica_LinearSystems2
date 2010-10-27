within Modelica_LinearSystems2.Controller;
block DAconverter "Digital to analog converter (including zero order hold)"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0)
    "Number of bits (=0 means no quantization error)";
  parameter Boolean unitDelay=true
    "= true, if one sample period delay, = false, if computing time not modelled";
  extends Interfaces.PartialSISO_equality;
protected
  Internal.DiscreteDAconverter discretePart(
    y_max=y_max,
    y_min=y_min,
    bits=bits,
    unitDelay=unitDelay,
    sampleFactor=sampleFactor) if
          not continuous "AD converter";
equation
  if continuous then
    y = if u > y_max then y_max else if u < y_min then y_min else u;
  else
    connect(u, discretePart.u);
    connect(y, discretePart.y);
  end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}), graphics={
        Text(
          extent={{30,-34},{96,-94}},
          lineColor={0,0,127},
          textString="A"),
        Line(points={{-100,-100},{100,100}}, color={0,0,127}),
        Text(
          extent={{-96,60},{-24,8}},
          lineColor={0,0,127},
          textString="D"),
        Line(
          points={{-26,24},{18,-20}},
          color={0,0,0},
          smooth=Smooth.None),
        Polygon(
          points={{14,-24},{22,-16},{26,-28},{14,-24}},
          lineColor={0,0,0},
          smooth=Smooth.None),
        Text(
          extent={{-92,92},{84,60}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}));
end DAconverter;
