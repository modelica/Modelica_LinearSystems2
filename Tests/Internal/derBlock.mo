within Modelica_LinearSystems2.Tests.Internal;
model derBlock
  annotation (uses(Modelica(version="3.1"), Modelica_LinearSystems2(version=
            "2.0")), Diagram(coordinateSystem(preserveAspectRatio=true, extent=
            {{-100,-100},{100,100}}), graphics));
  Modelica.Blocks.Sources.Sine sine(freqHz=0.1) 
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica_LinearSystems2.Controller.Sampler sampler 
    annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
  Modelica_LinearSystems2.Controller.Derivative derivative(T=0.1) 
    annotation (Placement(transformation(extent={{0,20},{20,40}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(
    blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Discrete,
    methodType=Modelica_LinearSystems2.Types.Method.ExplicitEuler,
    sampleTime=0.1) 
    annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
equation
  connect(sampler.u, sine.y) annotation (Line(
      points={{-42,30},{-59,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(derivative.u, sampler.y) annotation (Line(
      points={{-2,30},{-19,30}},
      color={0,0,127},
      smooth=Smooth.None));
end derBlock;
