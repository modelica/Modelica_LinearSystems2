within Modelica_LinearSystems2.Controllers.Interfaces;
partial block PartialSISO2
  "Partial Single Input Single Output (continuous or discrete) control block of Controllers library"
  extends PartialSampledBlock;
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block"
    annotation(Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block"
    annotation(Placement(transformation(extent={{100,-10},{120,10}})));
protected
  Internal.DiscreteStateSpace2 discretePart(
    methodType=methodType,
    sampleFactor=sampleFactor,
    init=init) if not continuous "Discretized SISO system";

equation
  connect(u, discretePart.u[1]);

end PartialSISO2;
