within Modelica_LinearSystems2.Controller;
block Sampler
  "Sample the input signal if discrete block or y=u if continuous block"
  extends Interfaces.PartialSISO_equality;

protected
  Internal.DiscreteSampler discretePart(sampleFactor=sampleFactor) if not continuous
    "Discrete sampler";
equation
  if continuous then
     y = u;
  else
    connect(u,discretePart.u);
    connect(y,discretePart.y);
  end if;
  annotation (
   Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Ellipse(
          extent={{45,-10},{25,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,0},{-45,0}}, color={0,0,127}),
        Line(points={{45,0},{100,0}}, color={0,0,127}),
        Line(points={{-35,0},{30,35}}, color={0,0,127}),
        Ellipse(
          extent={{-25,-10},{-45,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-90,-60},{90,-90}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-100,0},{-60,0}}),
        Line(points={{60,0},{100,0}}),
        Ellipse(
          extent={{45,-10},{25,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,0},{-45,0}}, color={0,0,127}),
        Line(points={{45,0},{100,0}}, color={0,0,127}),
        Line(points={{-35,0},{30,35}}, color={0,0,127}),
        Ellipse(
          extent={{-25,-10},{-45,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<html>
<p>
If <strong>discrete</strong> block, the output y is sampled according to sample time
sampleClock.sampleTime * sampleFactor, where sampleClock.sampleTime
is defined globally in the outer component sampleClock and
sampleFactor is an Integer parameter of component Sampler.
</p>
<p>
If <strong>continuous</strong> block, the output y is identical to the input u.
</p>
</html>"));
end Sampler;
