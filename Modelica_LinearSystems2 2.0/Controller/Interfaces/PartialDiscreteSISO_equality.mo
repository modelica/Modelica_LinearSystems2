within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialDiscreteSISO_equality
  "Common discrete part of blocks extending from PartialSISO_equality"
  extends Interfaces.PartialBlockIcon;

  parameter Integer sampleFactor(min=1)=1
    "Ts=sampleClock.sampleTime*sampleFactor";
  final parameter Modelica.SIunits.Time Ts = sampleClock.sampleTime*sampleFactor
    "Sample time" annotation(Hide=false);
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
    annotation (extent=[-140, -20; -100, 20]);
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
    annotation (extent=[100, -10; 120, 10]);

  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon,
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram,
    Documentation(info="<HTML>
</HTML>
"));
protected
  outer SampleClock sampleClock "Global options"                      annotation(Hide=true);

 // Derived quantities
  discrete Real u_sampled "Sampled continuous input signal u";
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(Hide=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(Hide=true);

equation
  if sampleClock.blockType == Types.BlockType.Continuous then
     // no sampling in sampleClock
     sampleTrigger = sample(Ts, Ts);
     ticks = 0;
  else
     when sampleClock.sampleTrigger then
        ticks = if pre(ticks) < sampleFactor then pre(ticks) + 1 else 1;
     end when;
     sampleTrigger = sampleClock.sampleTrigger and ticks >= sampleFactor;
  end if;

initial equation
  pre(ticks) = 0;
end PartialDiscreteSISO_equality;
