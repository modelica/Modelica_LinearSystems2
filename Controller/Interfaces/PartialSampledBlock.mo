within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialSampledBlock
  "Partial block of Sampled library (icon + default parameters)"

  extends Icons.PartialBlockIcon(cont=continuous);

  parameter Types.BlockTypeWithGlobalDefault blockType=Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block (Continuous/Discrete)" 
    annotation(Evaluate=true, Hide=true,Dialog(tab="Advanced options"));
  final parameter Boolean continuous = blockType == Types.BlockTypeWithGlobalDefault.Continuous or 
                                 blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and 
                                 sampleClock.blockType == Types.BlockType.Continuous
    "= true, if continuous block, otherwise discrete block";
  parameter Types.MethodWithGlobalDefault methodType=Types.MethodWithGlobalDefault.UseSampleClockOption if 
       not continuous "Type of discretization if discrete block" 
     annotation(Evaluate=true, Hide=true,Dialog(tab="Advanced options",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  final parameter Types.Init init=if initType == Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.UseSampleClockOption then 
            sampleClock.initType else initType
    "Type of initialization (no init/InitialState/SteadyState)"    annotation(Evaluate=true);

  parameter Integer sampleFactor(min=1)=1 if  not continuous
    "Ts=sampleClock.sampleTime*sampleFactor" 
     annotation(Dialog(tab="Advanced options",
                enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));

  parameter Types.InitWithGlobalDefault initType=Types.InitWithGlobalDefault.UseSampleClockOption
    "Type of initialization (no init/initial/steady state/output)" 
    annotation(Evaluate=true, Hide=true,  Dialog(tab="Advanced options"));

protected
  outer SampleClock sampleClock "Global options";
end PartialSampledBlock;
