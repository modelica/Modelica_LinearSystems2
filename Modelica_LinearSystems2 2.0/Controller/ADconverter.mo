within Modelica_LinearSystems2.Controller;
block ADconverter "Analog to digital converter (including sampler)"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0)
    "Number of bits (=0 means no quantization error)";
  extends Interfaces.PartialSISO_equality;

   annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-100,-100},{100,100}}, color={0,0,127}),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="bits=%bits"),
        Text(
          extent={{-98,98},{-18,18}},
          lineColor={0,0,127},
          textString="A"),
        Text(
          extent={{18,-18},{98,-98}},
          lineColor={0,0,127},
          textString="D"),
        Line(
          points={{-28,20},{16,-24}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{12,-28},{20,-20},{24,-32},{12,-28}},
          lineColor={0,0,255},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p>
If <b>discrete</b> block, the output y is sampled according to sample time
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
If <b>continuous</b> block, the output y is identical to the input u,
but is limited by y_min and y_max.
</p>
</html>"),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Window(
      x=0.45,
      y=0.01,
      width=0.35,
      height=0.49));

protected
  Internal.DiscreteADconverter discretePart(y_max=y_max, y_min=y_min, bits=bits,
      sampleFactor=sampleFactor) if
          not continuous "AD converter";

equation
   if continuous then
      y = if u > y_max then y_max else if u < y_min then y_min else u;
else
    connect(u, discretePart.u);
    connect(y, discretePart.y);
   end if;
end ADconverter;
