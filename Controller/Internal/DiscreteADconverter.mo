within Modelica_LinearSystems2.Controller.Internal;
block DiscreteADconverter "AD converter as discrete block"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0)
    "Number of bits (=0 means no quantization error)";
  extends Interfaces.PartialDiscreteSISO_equality;

protected
  parameter Real quantization=if bits > 0 then ((y_max - y_min)/2^bits) else 0;
  Real y_bound "Bounded output"
                               annotation(HideResult=true);
  discrete Real y_sampled "Sampled output" annotation(HideResult=true);
equation
  when {initial(), sampleTrigger} then
     u_sampled = u;
     y_bound = if u > y_max then y_max else if u < y_min then y_min else u;
     y_sampled = if bits > 0 then quantization*floor(abs(y_bound/quantization) + 0.5)*(
                if y_bound >= 0 then +1 else -1) else y_bound;
  end when;
  y = y_sampled;
  annotation (
    Documentation(info="<html>
</html>"), Icon(graphics={Line(
          points={{-100,-100},{100,100}},
          color={95,95,95},
          smooth=Smooth.None), Text(
          extent={{-90,90},{-10,10}},
          lineColor={95,95,95},
          textString="A"),     Text(
          extent={{10,-10},{90,-90}},
          lineColor={95,95,95},
          textString="D")}));
end DiscreteADconverter;
