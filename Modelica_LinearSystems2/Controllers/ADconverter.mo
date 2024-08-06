within Modelica_LinearSystems2.Controllers;
block ADconverter "Analog to digital converter (including sampler)"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0) "Number of bits (=0 : no quantization error)";
  extends Interfaces.PartialSISO_equality;

protected
  Internal.DiscreteADconverter discretePart(
    y_max=y_max,
    y_min=y_min,
    bits=bits,
    sampleFactor=sampleFactor) if not continuous "AD converter";

equation
  if continuous then
    y = if u > y_max then y_max else if u < y_min then y_min else u;
  else
    connect(u, discretePart.u);
    connect(y, discretePart.y);
  end if;
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Line(points={{-100,-100},{100,100}}, color={0,0,127}),
        Text(
          extent={{-94,60},{-30,20}},
          textColor={0,0,127},
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
          extent={{26,-10},{90,-50}},
          textColor={0,0,127},
          textString="D"),
        Text(
          extent={{-100,94},{100,64}},
          textColor={0,0,0},
          textString="%sampleFactor")}),
    Documentation(
      info="<html>
<p>
If <strong>discrete</strong> block, the output y is sampled according to sample time
sampleClock.sampleTime * sampleFactor, where sampleClock.sampleTime
is defined globally in the outer component sampleClock and
sampleFactor is an Integer parameter of component Sampler.
</p>
<p>
The sampled output signal is computed by limiting the input u with the
provided y_min and y_max borders and by rounding according to the
provided precision of the AD converter defined via parameter bits
(e.g. bits = 12 is the precision of simple AD converters).
</p>
<p>
If <strong>continuous</strong> block, the output y is identical to the input u,
but is limited by y_min and y_max.
</p>
</html>"));
end ADconverter;
