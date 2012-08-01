within Modelica_LinearSystems2.Controller.Internal;
block DiscreteSampler "Sample the input signal"
  extends Interfaces.PartialDiscreteSISO_equality;

equation
  when {initial(), sampleTrigger} then
      u_sampled = u;
  end when;

  y = u_sampled;
  annotation (
    Documentation(info="<html>
</html>"), Icon(graphics={
        Line(
          points={{-100,0},{-36,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{36,0},{100,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{-36,0},{28,36}},
          color={0,0,127},
          smooth=Smooth.None),
                   Ellipse(
          extent={{-26, -10},{-46, 10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{46,-10},{26,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}),
    Diagram(graphics={
        Line(
          points={{-100,0},{-36,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{36,0},{100,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{-36,0},{28,36}},
          color={0,0,127},
          smooth=Smooth.None),
                   Ellipse(
          extent={{-26,-10},{-46,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{46,-10},{26,10}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end DiscreteSampler;
