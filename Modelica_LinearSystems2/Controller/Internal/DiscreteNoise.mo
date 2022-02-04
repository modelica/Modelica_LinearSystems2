within Modelica_LinearSystems2.Controller.Internal;
block DiscreteNoise "Uniform distributed noise for discrete blocks"
  extends Controller.Icons.PartialBlockIcon(cont=false);

  parameter Real y_min "Lower limit of noise band";
  parameter Real y_max "Upper limit of noise band";
  parameter Integer firstSeed[3](each min=0, each max=255) = {23,87,187}
    "Integer[3] defining random sequence; required element range: 0..255";
  parameter Integer sampleFactor(min=1)=1
    "Noise sample time = sampleClock.sampleTime*sampleFactor";
  final parameter Modelica.Units.SI.Time Ts = sampleClock.sampleTime*sampleFactor
    "Sample time";
  Modelica.Blocks.Interfaces.RealOutput y
    "Noise output signal in the range [y_min .. y_max]"
   annotation (Placement(transformation(extent={{100,-10},{120,10}})));

protected
  outer SampleClock sampleClock "Global options";
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(HideResult=true);
  Integer seedState[3] "State of seed"
    annotation(HideResult=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(HideResult=true);
  discrete Real noise "Noise in the range 0..1"
    annotation(HideResult=true);
  discrete Real y_sampled "Sampled output" annotation(HideResult=true);
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
    Documentation(info="<html>
</html>"),
    Icon(
      graphics={
        Line(
          points={{-84,76},{-84,-92}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{86,-84},{-94,-84}},
          color={175,175,175},
          smooth=Smooth.None),
        Polygon(
          points={{-84,90},{-92,68},{-76,68},{-84,90}},
          lineColor={175,175,175},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{90,-84},{68,-92},{68,-76},{90,-84}},
          lineColor={175,175,175},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-35,13},{-35,-47},{-25,-47},{-25,-29},{-15,-29},{-15,-57},{-5,-57},
          {-5,25},{1,25},{1,39},{7,39},{7,-17},{17,-17},{17,-5},{23,-5},{23,-35},
          { 33,-35},{33,37},{43,37},{43,3},{51,3},{51,-63},{61,-63}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{-81,-29},{-67,-29},{-67,-13},{-59,-13},{-59,-61},{-51,-61},{-51,-39},
          {-43,-39},{-43,45},{-35,45},{-35,13}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{-90,-70},{84,-70}},
          color={255,0,0},
          smooth=Smooth.None),
        Line(
          points={{-89,50},{85,50}},
          color={255,0,0},
          smooth=Smooth.None),
        Text(
          extent={{-50,90},{50,60}},
          lineColor={95,95,95},
          textString="noise"),
        Text(
          extent={{-130,-100},{130,-130}},
          lineColor={0,0,0},
          textString="[%y_min .. %y_max]")}));
end DiscreteNoise;
