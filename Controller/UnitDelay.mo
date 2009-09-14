within Modelica_LinearSystems2.Controller;
block UnitDelay
  "Delay the input by a multiple of the base sample time if discrete block or y=u if continuous block"
  extends Interfaces.PartialSISO_equality;

  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Text(
          extent={{-94,10},{86,90}},
          lineColor={0,0,127},
          textString="1"),
        Line(points={{-66,0},{72,0}}, color={0,0,127}),
        Text(
          extent={{-94,-10},{86,-90}},
          lineColor={0,0,127},
          textString="z")}),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
<p>
If <b>discrete</b> block, the output y is sampled and is the value
of the sampled input signal u at the previous sample instant, where
sample time = sampleClock.sampleTime * sampleFactor and 
sampleClock.sampleTime is defined globally in the outer component 
sampleClock and sampleFactor is an Integer parameter of component UnitDelay.
</p>
<p>
If <b>continuous</b> block, the output y is identical to the input u.
</p>
</html>"));
protected
  Internal.DiscreteUnitDelay discretePart(sampleFactor=sampleFactor) if 
          not continuous "Discrete unit delay";
equation
  if continuous then
     y = u;

  end if;
   connect(u, discretePart.u);
    connect(y, discretePart.y);

end UnitDelay;
