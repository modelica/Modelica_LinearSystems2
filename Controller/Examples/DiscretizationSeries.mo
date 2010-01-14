within Modelica_LinearSystems2.Controller.Examples;
model DiscretizationSeries
  "Demonstrates the discretization methods for a series connection"
  extends Modelica.Icons.Example;
  parameter Modelica.SIunits.Time T1=0.2;
  parameter Modelica.SIunits.Time T2=0.15;
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(sampleTime=
        0.05) 
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  Modelica.Blocks.Sources.Step step(
    height=1.2,
    offset=0.2,
    startTime=0.1)                   annotation (extent=[-80,40; -60,60],
      Placement(transformation(extent={{-80,20},{-60,40}})));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}), graphics), Documentation(info="<html>
<p>
Demonstrates the different discretization methods by simulating the step
response of a second order system as continuous system and as discrete system
with the supported discretization methods. The step starts with an offset at 0.1 s
to demonstrate the steady-state initialization.
</p>

</html>"),
    experiment(StopTime=1.5, Tolerance=1e-006),
    experimentSetupOutput);

  FirstOrder S1(blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous,
      T=T1) annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
  FirstOrder S2(blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous,
      T=T2) annotation (Placement(transformation(extent={{0,20},{20,40}})));
  ZerosAndPoles S12(system=
        Modelica_LinearSystems2.ZerosAndPoles.'constructor'.fromZerosAndPoles(p=
         {Modelica_LinearSystems2.Math.Complex(re=-1/T1, im=0),
        Modelica_LinearSystems2.Math.Complex(re=-1/T2, im=0)}, k=1/(T1*T2))) 
    annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
  Modelica.Blocks.Math.Feedback diff 
    annotation (Placement(transformation(extent={{30,20},{50,40}})));
equation
  connect(step.y, S1.u)            annotation (Line(
      points={{-59,30},{-42,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S1.y, S2.u)                       annotation (Line(
      points={{-19,30},{-2,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, S12.u) annotation (Line(
      points={{-59,30},{-50,30},{-50,0},{-22,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S2.y, diff.u1) annotation (Line(
      points={{21,30},{32,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(S12.y, diff.u2) annotation (Line(
      points={{1,0},{40,0},{40,22}},
      color={0,0,127},
      smooth=Smooth.None));
end DiscretizationSeries;
