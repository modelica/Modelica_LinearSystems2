within Modelica_LinearSystems2.Controller.Templates;
partial model SimpleStateSpaceControl
  "Template for simple state feedback controllers with an optional pre-filter"

  MatrixGain feedbackMatrix 
    annotation (Placement(transformation(extent={{40,-50},{20,-30}})));
  MatrixGain preFilter(K=[0]) 
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  Modelica.Blocks.Math.Feedback feedback[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
  replaceable Modelica_LinearSystems2.Controller.Templates.PartialPlantMIMO
    plant(n=feedbackMatrix.nin, m=feedbackMatrix.nout) constrainedby
    Modelica_LinearSystems2.Controller.Templates.PartialPlantMIMO 
    annotation (Placement(transformation(extent={{60,-10},{80,10}})));
  Sampler samplerPreFilter[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{-35,-5},{-25,5}})));
  Sampler samplerFeedback[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{5,-45},{-5,-35}})));
  Sampler samplerOut[feedbackMatrix.nin] 
    annotation (Placement(transformation(extent={{65,-45},{55,-35}})));
  inner SampleClock sampleClock 
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  connect(feedback.y, plant.u) annotation (Line(
      points={{-1,0},{58,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedbackMatrix.u, samplerOut.y) annotation (Line(
      points={{42,-40},{54.5,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.u, feedbackMatrix.y) annotation (Line(
      points={{6,-40},{19,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.y, feedback.u2) annotation (Line(
      points={{-5.5,-40},{-10,-40},{-10,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preFilter.y, samplerPreFilter.u) annotation (Line(
      points={{-39,0},{-36,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerPreFilter.y, feedback.u1) annotation (Line(
      points={{-24.5,0},{-18,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, samplerOut.u) annotation (Line(
      points={{70,-11},{70,-40},{66,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (defaultComponentName="controller",Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-140,
            -100},{140,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-140,-100},{140,100}})));
end SimpleStateSpaceControl;
