within Modelica_LinearSystems2.Controller.Types;
type BlockTypeWithGlobalDefault = enumeration(
    Continuous "Continuous block",
    Discrete "Discrete block",
    UseSampleClockOption "Use blockType defined in sampleClock component")
  "Enumeration of block types (continuous, discrete, sampleBlock.blocktype)"
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));
