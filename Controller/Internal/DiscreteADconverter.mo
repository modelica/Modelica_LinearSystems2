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
                               annotation(Hide=true);
  discrete Real y_sampled "Sampled output" annotation(Hide=true);
equation
  when {initial(), sampleTrigger} then
     u_sampled = u;
     y_bound = if u > y_max then y_max else if u < y_min then y_min else u;
     y_sampled = if bits > 0 then quantization*floor(abs(y_bound/quantization) + 0.5)*(
                if y_bound >= 0 then +1 else -1) else y_bound;
  end when;
  y = y_sampled;
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(
      Text(
        extent=[-90,90; -10,10],
        string="A",
          style(color=74, rgbcolor={0,0,127})),
      Line(points=[-100,-100; 100,100], style(color=74, rgbcolor={0,0,127})),
      Text(
        extent=[10,-10; 90,-90],
        string="D",
          style(color=74, rgbcolor={0,0,127}))),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram,
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteADconverter;
