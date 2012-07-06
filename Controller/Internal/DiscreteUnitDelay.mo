within Modelica_LinearSystems2.Controller.Internal;
block DiscreteUnitDelay "Delay the input signal by one sample instant"
  extends Interfaces.PartialDiscreteSISO_equality;

protected
  discrete Real y_sampled "Sampled output" annotation(Hide=true);
equation
  when {initial(), sampleTrigger} then
     u_sampled = u;
     y_sampled = pre(u_sampled);
  end when;
  y = y_sampled;
  annotation (
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteUnitDelay;
