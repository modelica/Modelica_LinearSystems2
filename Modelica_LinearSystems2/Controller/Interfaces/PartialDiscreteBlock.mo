within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialDiscreteBlock
  "Partial discrete block, i.e. no continuous alternative"
  extends Controller.Icons.PartialBlockIcon(cont=false);

  parameter Types.InitWithGlobalDefault initType=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (NoInit/SteadyState/InitialState/InitialOutput)"
    annotation(Evaluate=true, HideResult=true,  Dialog(tab="Advanced options"));
  final parameter Types.Init init=Modelica_LinearSystems2.Controller.Internal.convertToInit(initType, sampleClock.initType)
    "Type of initialization (NoInit/SteadyState/InitialState/InitialOutput)" annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1) = 1
    "Sample factor (Ts = sampleFactor * sampleClock.sampleTime)"
    annotation(HideResult=true);
  final parameter Modelica.Units.SI.Time Ts=sampleClock.sampleTime*sampleFactor
    "Sample time" annotation (HideResult=false);
protected
  Integer ticks
    "Actual number of base samples starting from the last sample time instant" annotation(HideResult=true);
    Boolean sampleTrigger "Triggers next sample time" annotation(HideResult=true);
  outer SampleClock sampleClock "Global options";

initial equation
  pre(ticks) = 0;
equation
  when sampleClock.sampleTrigger then
    ticks = if pre(ticks) < sampleFactor then pre(ticks) + 1 else 1;
  end when;
  sampleTrigger = sampleClock.sampleTrigger and ticks >= sampleFactor;
end PartialDiscreteBlock;
