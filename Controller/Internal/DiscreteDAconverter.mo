within Modelica_LinearSystems2.Controller.Internal;
block DiscreteDAconverter
  "Digital to analog converter as discrete block (including zero order hold)"
  parameter Real y_max "Upper limit of output signal";
  parameter Real y_min "Lower limit of output signal";
  parameter Integer bits(min=0)
    "Number of bits (=0 means no quantization error)";
  parameter Boolean unitDelay = true
    "= true, if one sample period delay, = false, if computing time not modelled";
  extends Interfaces.PartialDiscreteSISO_equality;
protected
  parameter Real quantization=if bits > 0 then ((y_max - y_min)/2^bits) else 0;
  discrete Real y_bound "Bounded output" 
                                annotation(Hide=true);
  discrete Real y_sampled "Sampled output" 
                                  annotation(Hide=true);
  discrete Real y_delaySampled
    "Sampled output with a delay of one sample period"                            annotation(Hide=true);
equation
  when {initial(),sampleTrigger} then
     u_sampled = u;
     y_bound = if u > y_max then y_max else if u < y_min then y_min else u;
     y_sampled = if bits > 0 then quantization*floor(abs(y_bound/quantization) + 0.5)
                 *(if y_bound >= 0 then +1 else -1) else y_bound;
     if unitDelay then
        y_delaySampled = if initial() then y_sampled else pre(y_sampled);
     else
        y_delaySampled = y_sampled;
     end if;
  end when;

  y = y_delaySampled;
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(
      Text(
        extent=[-90, 90; -10, 10],
        string="D",
        style(color=74, rgbcolor={0,0,127})),
      Text(
        extent=[10, -10; 90, -90],
        string="A",
        style(color=74, rgbcolor={0,0,127})),
      Line(points=[-100, -100; 100, 100], style(color=74, rgbcolor={0,0,127})),
      Text(
        extent=[-150, -150; 150, -110],
        string="bits=%bits",
        style(color=0))),
    Diagram(
      Text(
        extent=[-90, 90; -10, 10],
        string="D",
        style(color=74, rgbcolor={0,0,127})),
      Text(
        extent=[10, -10; 90, -90],
        string="A",
        style(color=74, rgbcolor={0,0,127})),
      Line(points=[-100, -100; 100, 100], style(color=74, rgbcolor={0,0,127})),
      Rectangle(extent=[-100, 100; 100, -100])),
    Documentation(info="<HTML>
</HTML>
"), Window(
      x=0.45,
      y=0.01,
      width=0.35,
      height=0.49));
end DiscreteDAconverter;
