within Modelica_LinearSystems2.Controllers.Interfaces;
partial block PartialSISO_equality
  "Partial Single Input Single Output (continuous or discrete) control block of Controllers library where usually y=u if continuous block"
  extends Icons.PartialBlockIcon(cont=continuous);
  import Modelica_LinearSystems2.Controllers.Types;

  parameter Types.BlockTypeWithGlobalDefault blockType=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block"
    annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(
        __Dymola_compact=true,
        __Dymola_descriptionLabel=true),
      choices(__Dymola_radioButtons=true, choice=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Continuous
        "Continuous",
        choice=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Discrete
        "Discrete",
        choice=Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.UseSampleClockOption
        "Dependent on sampleClock"));
  final parameter Boolean continuous = blockType == Types.BlockTypeWithGlobalDefault.Continuous or
                                 blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Continuous
    "True, if continuous block, otherwise discrete block";
  parameter Integer sampleFactor(min=1)=1
    "Sample factor for sample time (Ts = sampleFactor * sampleClock.sampleTime)"
     annotation (Dialog(enable=blockType<>Modelica_LinearSystems2.Controllers.Types.BlockTypeWithGlobalDefault.Continuous,
     group="Discrete block parameters"));
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
     annotation(Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
     annotation(Placement(transformation(extent={{100,-10},{120,10}})));

protected
  outer SampleClock sampleClock "Global options";
  annotation (
    Documentation(info="<html>
<p>
If <b>discrete</b> block, the output y is sampled according to sample time
sampleClock.sampleTime * sampleFactor, where sampleClock.sampleTime
is defined globally in the outer component sampleClock and
sampleFactor is an Integer parameter of component Sampler.
</p>
<p>
If <b>continuous</b> block, the output y is identical to the input u.
</p>
</html>"));
end PartialSISO_equality;
