within Modelica_LinearSystems2.Controller.Internal;
block DiscreteNoise "Uniform distributed noise for discrete blocks"
  extends Controller.Icons.PartialBlockIcon(cont=false);

  parameter Real y_min "Lower limit of noise band";
  parameter Real y_max "Upper limit of noise band";
  parameter Integer firstSeed[3](each min=0, each max=255) = {23,87,187}
    "Integer[3] defining random sequence; required element range: 0..255";
  parameter Integer sampleFactor(min=1)=1
    "Noise sample time = sampleClock.sampleTime*sampleFactor";
  final parameter Modelica.SIunits.Time Ts = sampleClock.sampleTime*sampleFactor
    "Sample time" annotation(Hide=false);
  Modelica.Blocks.Interfaces.RealOutput y
    "Noise output signal in the range [y_min .. y_max]"
    annotation (extent=[100, -10; 120, 10]);

protected
  outer SampleClock sampleClock "Global options";
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(Hide=true);
  Integer seedState[3] "State of seed"
                       annotation(Hide=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(Hide=true);
  discrete Real noise "Noise in the range 0..1"
                                       annotation(Hide=true);
  discrete Real y_sampled "Sampled output" annotation(Hide=true);
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

  when {initial(), sampleTrigger} then
     (noise,seedState) = random(pre(seedState));
      y_sampled = y_min + (y_max - y_min)*noise;
  end when;
  y = y_sampled;

initial equation
  pre(ticks) = 0;
  pre(seedState) = firstSeed;
  annotation (
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteNoise;
