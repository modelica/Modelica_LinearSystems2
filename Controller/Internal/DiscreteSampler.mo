within Modelica_LinearSystems2.Controller.Internal;
block DiscreteSampler "Sample the input signal"
  extends Interfaces.PartialDiscreteSISO_equality;

equation
  when {initial(), sampleTrigger} then
      u_sampled = u;
  end when;

  y = u_sampled;
  annotation (
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteSampler;
