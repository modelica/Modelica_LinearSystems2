within Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities;
model ControlledSISO1 "Controlled SISO system"
  extends Modelica.Blocks.Interfaces.SISO;
  parameter Real k=1 "Input signal gain";

  Controllers.ZerosAndPoles zerosAndPoles(
    system(
      n2=[2,2; 4,8],
      d2=[20,101; 22,122],
      n1=fill(0.0, 0),
      d1=fill(0.0, 0)))
    annotation (Placement(transformation(extent={{10,-10},{30,10}})));
  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Math.Gain gain(k=k)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  inner Controllers.SampleClock sampleClock
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
equation
  connect(zerosAndPoles.u, gain.y) annotation (Line(
      points={{8,0},{-9,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, feedback.y) annotation (Line(
      points={{-32,0},{-61,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.y, feedback.u2) annotation (Line(
      points={{31,0},{60,0},{60,-40},{-70,-40},{-70,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.u1, u) annotation (Line(
      points={{-78,0},{-92,0},{-92,0},{-120,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zerosAndPoles.y, y) annotation (Line(
      points={{31,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics), Documentation(info="<html>
<p>
Utility model in order to demonstrate the plotting of a root locus:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/ControlledSISO1a.png\"/></p>

<p>
This model consists of a linear, time invariant single-input, single-output plant &quot;zerosAndPoles&quot;
that is controlled by a P-Controller with a constant gain k. The pole/zero pattern of the plant
is shown in the next figure:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/ControlledSISO1b.png\"/></p>
</html>"));
end ControlledSISO1;
