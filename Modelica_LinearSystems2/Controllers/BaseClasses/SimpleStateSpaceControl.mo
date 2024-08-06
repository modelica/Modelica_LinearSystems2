within Modelica_LinearSystems2.Controllers.BaseClasses;
partial model SimpleStateSpaceControl
  "Template for a simple state feedback controller with an optional pre-filter"

  MatrixGain feedbackMatrix
    annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
  MatrixGain preFilter(K=[0])
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  Modelica.Blocks.Math.Feedback feedback[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{0,-10},{20,10}})));
  replaceable Modelica_LinearSystems2.Controllers.BaseClasses.PartialPlantMIMO plant(
    n=feedbackMatrix.nin,
    m=feedbackMatrix.nout) constrainedby Modelica_LinearSystems2.Controllers.BaseClasses.PartialPlantMIMO
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  Sampler samplerPreFilter[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{-15,-5},{-5,5}})));
  Sampler samplerFeedback[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{25,-45},{15,-35}})));
  Sampler samplerOut[feedbackMatrix.nin]
    annotation (Placement(transformation(extent={{85,-45},{75,-35}})));
  inner SampleClock sampleClock
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  connect(feedback.y, plant.u) annotation (Line(
      points={{19,0},{78,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedbackMatrix.u, samplerOut.y) annotation (Line(
      points={{62,-40},{74.5,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.u, feedbackMatrix.y) annotation (Line(
      points={{26,-40},{39,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.y, feedback.u2) annotation (Line(
      points={{14.5,-40},{10,-40},{10,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preFilter.y, samplerPreFilter.u) annotation (Line(
      points={{-19,0},{-16,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerPreFilter.y, feedback.u1) annotation (Line(
      points={{-4.5,0},{2,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, samplerOut.u) annotation (Line(
      points={{90,-11},{90,-40},{86,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (defaultComponentName="controller");
end SimpleStateSpaceControl;
