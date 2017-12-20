within Modelica_LinearSystems2.Utilities.Plants;
block SimpleLinearSystem "Simple linear test system with zeros and poles"
  extends Modelica.Blocks.Interfaces.SISO;

  Modelica.Blocks.Continuous.Filter filter(
    filterType=Modelica.Blocks.Types.FilterType.HighPass,
    order=3,
    f_cut=5,
    analogFilter=Modelica.Blocks.Types.AnalogFilter.Bessel)
             annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
equation
  connect(filter.u, u) annotation (Line(
      points={{-12,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter.y, y) annotation (Line(
      points={{11,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Documentation(info="<html>
<p>A simple linear SISO model of test system with zeros and poles.</p>
<!--placeholder-->
</html>"));
end SimpleLinearSystem;
