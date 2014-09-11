within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialSampledBlock
  "Partial block of Sampled library (icon + default parameters)"
  extends Icons.PartialBlockIcon(cont=continuous);
  import Modelica_LinearSystems2.Controller.Types;

  parameter Types.BlockTypeWithGlobalDefault blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block"
    annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(
        tab="Advanced options",
        __Dymola_compact=true,
        __Dymola_descriptionLabel=true),
      choices(__Dymola_radioButtons=true, choice=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous
        "Continuous",
        choice=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete
        "Discrete",
        choice=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption
        "Dependent on sampleClock"));
  final parameter Boolean continuous = blockType == Types.BlockTypeWithGlobalDefault.Continuous or
                                 blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Continuous
    "True, if continuous block, otherwise discrete block";
  parameter Types.MethodWithGlobalDefault methodType=Types.MethodWithGlobalDefault.UseSampleClockOption
    "Type of discretization if discrete block"
     annotation(Evaluate=true, HideResult=true,Dialog(tab="Advanced options",group="Discrete block parameters",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  final parameter Types.Init init=if initType == Types.InitWithGlobalDefault.UseSampleClockOption then
            sampleClock.initType else
            (if initType==Types.InitWithGlobalDefault.NoInit then Types.Init.NoInit else
            (if initType==Types.InitWithGlobalDefault.SteadyState then Types.Init.SteadyState else
            (if initType==Types.InitWithGlobalDefault.InitialState then Types.Init.InitialState else
            Types.Init.InitialOutput)))
    "Type of initialization (no init/InitialState/SteadyState)" annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1)=1
    "Sample factor (Ts = sampleFactor * sampleClock.sampleTime)"
     annotation(Dialog(tab="Advanced options",group="Discrete block parameters",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  parameter Types.InitWithGlobalDefault initType=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (no init/initial/steady state/output)"
    annotation(Evaluate=true, HideResult=true,  Dialog(tab="Advanced options"));

protected
  outer SampleClock sampleClock "Global options";
end PartialSampledBlock;
