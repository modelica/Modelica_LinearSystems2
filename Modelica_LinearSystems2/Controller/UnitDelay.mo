within Modelica_LinearSystems2.Controller;
block UnitDelay
  "Delay the input by a multiple of the base sample time if discrete block or y=u if continuous block"
  extends Interfaces.PartialSISO_equality;

protected
  Internal.DiscreteUnitDelay discretePart(sampleFactor=sampleFactor) if
          not continuous "Discrete unit delay";
equation
  if continuous then
    y = u;
  end if;

  connect(u, discretePart.u);
  connect(y, discretePart.y);

  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Text(
          extent={{-92,26},{88,92}},
          lineColor={0,0,127},
          textString="1"),
        Line(points={{-70,20},{68,20}},
          color={0,0,127}),
        Text(
          extent={{-86,22},{88,-62}},
          lineColor={0,0,127},
          textString="z"),
        Text(
          extent={{-90,-60},{90,-90}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Documentation(info="<html>
<p>
If <strong>discrete</strong> block, the output y is sampled and is the value
of the sampled input signal u at the previous sample instant, where
sample time = sampleClock.sampleTime * sampleFactor and
sampleClock.sampleTime is defined globally in the outer component
sampleClock and sampleFactor is an Integer parameter of component UnitDelay.
</p>
<p>
If <strong>continuous</strong> block, the output y is identical to the input u.
</p>
</html>"));
end UnitDelay;
