within Modelica_LinearSystems2.Controllers;
block DAconverter "Digital to analog converter (including zero order hold)"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0) "Number of bits (=0 : no quantization error)";
  parameter Boolean unitDelay=true
    "True, if one sample period delay, otherwise computing time not modelled";
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
        Line(points={{-100,-100},{100,100}}, color={0,0,127}),
        Text(
          extent={{-94,60},{-30,20}},
          lineColor={0,0,127},
          textString="D"),
        Text(
          extent={{26,-10},{90,-50}},
          lineColor={0,0,127},
          textString="A"),
        Line(
          points={{-28,28},{24,-24}},
          color={0,0,0},
          smooth=Smooth.None),
        Polygon(
          points={{12,-20},{20,-12},{24,-24},{12,-20}},
          lineColor={0,0,0},
          smooth=Smooth.None),
        Text(
          extent={{-100,90},{80,60}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}));
end DAconverter;
