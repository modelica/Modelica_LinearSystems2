within Modelica_LinearSystems2.Controller.Internal;
block DiscreteUnitDelay "Delay the input signal by one sample instant"
  extends Interfaces.PartialDiscreteSISO_equality;

protected
  discrete Real y_sampled "Sampled output" annotation(HideResult=true);
initial equation
  pre(u_sampled) = 0.0;
equation
  when {initial(), sampleTrigger} then
     u_sampled = u;
     y_sampled = pre(u_sampled);
  end when;
  y = y_sampled;
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
          extent={{-26,-10},{-46,10}},
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
end DiscreteUnitDelay;
