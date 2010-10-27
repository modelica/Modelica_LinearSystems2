within Modelica_LinearSystems2.Controller;
block Sampler
  "Sample the input signal if discrete block or y=u if continuous block"
  extends Interfaces.PartialSISO_equality;

protected
  Internal.DiscreteSampler discretePart(sampleFactor=sampleFactor) if
          not continuous "Discrete sampler";
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
          extent={{-25,-10},{-45,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{45,-10},{25,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,0},{-45,0}}, color={0,0,127}),
        Line(points={{45,0},{100,0}}, color={0,0,127}),
        Line(points={{-35,0},{30,35}}, color={0,0,127}),
        Text(
          extent={{-98,-30},{98,-64}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-100,0},{-60,0}}),
        Line(points={{60,0},{100,0}}),
        Ellipse(
          extent={{-25,-10},{-45,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{45,-10},{25,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,0},{-45,0}}, color={0,0,127}),
        Line(points={{45,0},{100,0}}, color={0,0,127}),
        Line(points={{-35,0},{30,35}}, color={0,0,127})}),
    Documentation(info="<html>
<p>
If <b>discrete</b> block, the output y is sampled according to sample time
sampleClock.sampleTime * sampleFactor, where sampleClock.sampleTime
is defined globally in the outer component sampleClock and
sampleFactor is an Integer parameter of component Sampler.
</p>
<p>
If <b>continuous</b> block, the output y is identical to the input u.
</p>
</html>"));
end Sampler;
