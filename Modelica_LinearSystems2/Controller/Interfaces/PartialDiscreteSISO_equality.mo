within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialDiscreteSISO_equality
  "Common discrete part of blocks extending from PartialSISO_equality"
  extends Controller.Icons.PartialBlockIcon(cont=false);

  parameter Integer sampleFactor(min=1)=1
    "Sample factor (Ts = sampleFactor * sampleClock.sampleTime)";
  final parameter Modelica.Units.SI.Time Ts = sampleClock.sampleTime*sampleFactor
    "Sample time" annotation (HideResult=false);
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
     annotation(Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
     annotation(Placement(transformation(extent={{100,-10},{120,10}})));

protected
  outer SampleClock sampleClock "Global options"                      annotation(HideResult=true);

 // Derived quantities
  discrete Real u_sampled "Sampled continuous input signal u";
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(HideResult=true);
  Boolean sampleTrigger "Triggers next sample time" annotation(HideResult=true);

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
  annotation (
    Documentation(info="<html>
</html>"));
end PartialDiscreteSISO_equality;
