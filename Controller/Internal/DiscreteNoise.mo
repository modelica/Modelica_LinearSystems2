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
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[1,1],
      component=[20, 20],
        scale=0),
    Icon(
      Polygon(points=[-81,90; -89,68; -73,68; -81,90],     style(
          color=8,
          fillColor=8,
          fillPattern=1)),
      Line(points=[-81,78; -81,-90],   style(color=8)),
      Line(points=[-90,-23; 82,-23],   style(color=8)),
      Polygon(points=[91,-22; 69,-14; 69,-30; 91,-22],     style(
          color=8,
          fillColor=8,
          fillPattern=1)),
      Text(
        extent=[-24,95; 94,56],
        style(color=8),
          string="noise"),
      Line(points=[-35,13; -35,-47; -25,-47; -25,-29; -15,-29; -15,-57; -5,
              -57; -5,25; 1,25; 1,39; 7,39; 7,-17; 17,-17; 17,-5; 23,-5; 23,
              -35; 33,-35; 33,37; 43,37; 43,3; 51,3; 51,-63; 61,-63], style(
              color=74, rgbcolor={0,0,127})),
      Line(points=[-81,-29; -67,-29; -67,-13; -59,-13; -59,-61; -51,-61; -51,
              -39; -43,-39; -43,45; -35,45; -35,13], style(color=74, rgbcolor=
               {0,0,127})),
        Line(points=[-90,-70; 84,-70], style(color=1, rgbcolor={255,0,0})),
        Line(points=[-89,50; 85,50], style(color=1, rgbcolor={255,0,0})),
        Text(
          extent=[-145,-105; 144,-133],
          string="[%y_min .. %y_max]",
          style(color=0, rgbcolor={0,0,0}))),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram,
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteNoise;
