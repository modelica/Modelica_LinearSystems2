within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialSISO_equality
  "Partial Single Input Single Output (continuous or discrete) control block of Controller library where usually y=u if continuous block"
  extends Icons.PartialBlockIcon(cont=continuous);
  import Modelica_LinearSystems2.Controller.Types;
  parameter Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault
    blockType =                                                                           Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block (Continuous/Discrete)"
    annotation(Evaluate=true, Hide=true);
  final parameter Boolean continuous = blockType == Types.BlockTypeWithGlobalDefault.Continuous or
                                 blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption and
                                 sampleClock.blockType == Types.BlockType.Continuous
    "= true, if continuous block, otherwise discrete block";
  parameter Integer sampleFactor(min=1)=1 if not continuous
    "Sample time = sampleFactor * sampleClock.sampleTime"
     annotation (Dialog(enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
    annotation (extent=[-140, -20; -100, 20]);
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
    annotation (extent=[100, -10; 120, 10]);

protected
  outer SampleClock sampleClock "Global options";
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon,
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram,
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
