within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialDiscreteBlock
  "Partial discrete block, i.e. no continuous alternative"
  extends Controller.Icons.PartialBlockIcon(cont=false);

  parameter Types.InitWithGlobalDefault initType=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (no init/initial/steady state/output)"
    annotation(Evaluate=true, HideResult=true,  Dialog(tab="Advanced options"));
  final parameter Types.Init init=if initType == Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.UseSampleClockOption then
            sampleClock.initType else initType
    "Type of initialization (no init/InitialState/SteadyState)" annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1) = 1
    "Sample factor (Ts = sampleFactor * sampleClock.sampleTime)"
     annotation(HideResult=true);
  final parameter Modelica.SIunits.Time Ts=sampleClock.sampleTime*sampleFactor
    "Sample time" annotation(HideResult=false);
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
