within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialSISO
  "Partial Single Input Single Output (continuous or discrete) control block of Controller library"
  extends PartialSampledBlock;
  Modelica.Blocks.Interfaces.RealInput u
    "Continuous or discrete input signal of block" 
    annotation (extent=[-140, -20; -100, 20]);
  Modelica.Blocks.Interfaces.RealOutput y
    "Continuous or discrete output signal of block" 
    annotation (extent=[100, -10; 120, 10]);

protected
  Internal.DiscreteStateSpace discretePart(
    methodType=methodType,
    sampleFactor=sampleFactor,
    init=init) if not continuous "Discretized SISO system";

equation
 connect(u, discretePart.u[1]);

  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Window(
      x=0.01,
      y=0.24,
      width=0.7,
      height=0.72),
    Documentation(info="<HTML>
</HTML>
"), Diagram);
end PartialSISO;
