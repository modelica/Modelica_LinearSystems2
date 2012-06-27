within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialSampledBlock
  "Partial block of Sampled library (icon + default parameters)"
  extends Icons.PartialBlockIcon(cont=continuous);

  parameter Types.BlockTypeWithGlobalDefault blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block"
    annotation (
      Evaluate=true,
      Hide=true,
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
  parameter Types.MethodWithGlobalDefault methodType=Types.MethodWithGlobalDefault.UseSampleClockOption if
       not continuous "Type of discretization if discrete block"
     annotation(Evaluate=true, Hide=true,Dialog(tab="Advanced options",group="Discrete block parameters",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  final parameter Types.Init init=if initType == Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.UseSampleClockOption then
            sampleClock.initType else initType
    "Type of initialization (no init/InitialState/SteadyState)" annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1)=1 if not continuous
    "Sample factor (Ts = sampleFactor * sampleClock.sampleTime)"
     annotation(Dialog(tab="Advanced options",group="Discrete block parameters",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  parameter Types.InitWithGlobalDefault initType=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (no init/initial/steady state/output)"
    annotation(Evaluate=true, Hide=true,  Dialog(tab="Advanced options"));

protected
  outer SampleClock sampleClock "Global options";
end PartialSampledBlock;
