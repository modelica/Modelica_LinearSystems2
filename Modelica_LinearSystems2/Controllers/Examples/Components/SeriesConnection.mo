within Modelica_LinearSystems2.Controllers.Examples.Components;
model SeriesConnection "Series connection of two PT1 blocks"
  extends Modelica.Blocks.Interfaces.BlockIcon;
  parameter Types.BlockType blockType=Modelica_LinearSystems2.Controllers.Types.BlockType.Continuous
    "Type of Sampled blocks (Continuous or Discrete)";
  parameter Modelica_LinearSystems2.Utilities.Types.Method methodType=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method for discrete blocks";
  parameter Modelica.SIunits.Time sampleTime=0.05
    "Base sample time for discrete blocks";
  parameter Modelica.SIunits.Time T1=0.2 "Time constant of first PT1 block";
  parameter Modelica.SIunits.Time T2=0.15 "Time constant of second PT1 block";

  inner Modelica_LinearSystems2.Controllers.SampleClock sampleClock(
    blockType=blockType,
    methodType=methodType,
    sampleTime=sampleTime)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica.Blocks.Sources.Step step(
    height=1.2,
    offset=0.2,
    startTime=0.1) annotation (
      Placement(transformation(extent={{-70,-10},{-50,10}})));

  FirstOrder S1(T=T1)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  FirstOrder S2(T=T2)
    annotation (Placement(transformation(extent={{10,-10},{30,10}})));
  ZerosAndPoles S12(
    system=Modelica_LinearSystems2.ZerosAndPoles.'constructor'.fromZerosAndPoles(
      p={Complex(re=-1/T1, im=0), Complex(re=-1/T2, im=0)}, k=1/(T1*T2)))
    annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
  Modelica.Blocks.Math.Feedback diff
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));

equation
  connect(step.y, S1.u)            annotation (Line(
      points={{-49,0},{-32,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S1.y, S2.u)                       annotation (Line(
      points={{-9,0},{8,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, S12.u) annotation (Line(
      points={{-49,0},{-40,0},{-40,-30},{-12,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S2.y, diff.u1) annotation (Line(
      points={{31,0},{42,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S12.y, diff.u2) annotation (Line(
      points={{11,-30},{50,-30},{50,-8}},
      color={0,0,127},
      smooth=Smooth.None));

  annotation (
    Icon(graphics={
        Rectangle(extent={{-66,26},{-20,-10}}, lineColor={0,0,127}),
        Rectangle(extent={{20,26},{72,-10}}, lineColor={0,0,127}),
        Line(
          points={{-20,8},{20,8}},
          color={0,0,127},
          smooth=Smooth.None),
        Polygon(
          points={{20,8},{6,12},{6,4},{20,8}},
          lineColor={0,0,127},
          smooth=Smooth.None,
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid)}), Documentation(info="<html>
<p>Two PT1 blocks connected in serie.</p>
</html>"));
end SeriesConnection;
