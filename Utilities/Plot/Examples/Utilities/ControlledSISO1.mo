within Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities;
model ControlledSISO1
  extends Modelica.Blocks.Interfaces.SISO;
  parameter Real k=1;

  Controller.ZerosAndPoles zerosAndPoles(system(
      n2=[2,2; 4,8],
      d2=[20,101; 22,122],
      n1=fill(0.0, 0),
      d1=fill(0.0, 0)))
    annotation (Placement(transformation(extent={{8,-10},{28,10}})));
  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Math.Gain gain(k=k)
    annotation (Placement(transformation(extent={{-36,-10},{-16,10}})));
  inner Controller.SampleClock sampleClock
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
equation
  connect(zerosAndPoles.u, gain.y) annotation (Line(
      points={{6,0},{-15,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, feedback.y) annotation (Line(
      points={{-38,0},{-61,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.y, feedback.u2) annotation (Line(
      points={{29,0},{52,0},{52,-38},{-70,-38},{-70,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, u) annotation (Line(
      points={{-78,0},{-92,0},{-92,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.y, y) annotation (Line(
      points={{29,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end ControlledSISO1;
