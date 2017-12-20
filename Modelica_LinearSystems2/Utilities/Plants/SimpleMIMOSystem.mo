within Modelica_LinearSystems2.Utilities.Plants;
block SimpleMIMOSystem
  "Simple MIMO test system with two inputs, two outputs, zeros and poles"
  extends Modelica.Blocks.Interfaces.MIMO(nin=2,nout=2);

  Modelica.Blocks.Continuous.Filter filter2(
    filterType=Modelica.Blocks.Types.FilterType.HighPass,
    order=3,
    analogFilter=Modelica.Blocks.Types.AnalogFilter.Butterworth,
    f_cut=2) annotation (Placement(transformation(extent={{-12,20},{8,40}})));
  Modelica.Blocks.Continuous.Filter filter1(
    filterType=Modelica.Blocks.Types.FilterType.HighPass,
    order=3,
    f_cut=2,
    analogFilter=Modelica.Blocks.Types.AnalogFilter.Bessel)
    annotation (Placement(transformation(extent={{-12,-40},{8,-20}})));
equation
  connect(filter1.u, u[1]) annotation (Line(
      points={{-14,-30},{-30,-30},{-30,-2},{-120,-2},{-120,-10},{-120,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter2.u, u[2]) annotation (Line(
      points={{-14,30},{-30,30},{-30,2},{-120,2},{-120,10},{-120,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter1.y, y[1]) annotation (Line(
      points={{9,-30},{30,-30},{30,-2},{110,-2},{110,-5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter2.y, y[2]) annotation (Line(
      points={{9,30},{30,30},{30,2},{110,2},{110,5}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Documentation(info="<html>
<p>A simple model of MIMO test system with two inputs, two outputs, zeros and poles.</p>
<!--placeholder-->
</html>
"));
end SimpleMIMOSystem;
