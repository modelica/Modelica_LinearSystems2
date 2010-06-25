within Modelica_LinearSystems2.Controller.Templates;
partial model SimpleObserverStateSpaceControl
  "Represents the structure of a simple state feedback controller with observer and optional pre filter"

  MatrixGain feedbackMatrix 
    annotation (Placement(transformation(extent={{20,-50},{0,-30}})));
  MatrixGain preFilter(K=[0]) 
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Math.Feedback feedback[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  replaceable Internal.PlantTemplate plant(n=
        feedbackMatrix.nin, m=feedbackMatrix.nout) constrainedby
    Internal.PlantTemplate 
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  Sampler samplerPreFilter[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{-53,-5},{-43,5}})));
  Sampler samplerFeedback[feedbackMatrix.nout] 
    annotation (Placement(transformation(extent={{-15,-45},{-25,-35}})));
  Sampler samplerOut[observer.nout] 
    annotation (Placement(transformation(extent={{85,-51},{75,-41}})));
  inner SampleClock sampleClock 
    annotation (Placement(transformation(extent={{120,80},{140,100}})));
  Modelica_LinearSystems2.Controller.Templates.Internal.ObserverTemplate
    observer 
    annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
equation
  connect(feedback.y, plant.u) annotation (Line(
      points={{-21,0},{82,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.u, feedbackMatrix.y) annotation (Line(
      points={{-14,-40},{-1,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.y, feedback.u2) annotation (Line(
      points={{-25.5,-40},{-30,-40},{-30,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preFilter.y, samplerPreFilter.u) annotation (Line(
      points={{-59,0},{-54,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerPreFilter.y, feedback.u1) annotation (Line(
      points={{-42.5,0},{-38,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerOut.y, observer.y)         annotation (Line(
      points={{74.5,-46},{62,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, samplerOut.u) annotation (Line(
      points={{90,-11},{90,-46},{86,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(observer.x_estimated, feedbackMatrix.u) annotation (Line(
      points={{39,-40},{22,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(observer.u, feedback.y) annotation (Line(
      points={{62,-34},{68,-34},{68,0},{-21,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-140,
            -100},{140,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-140,-100},{140,100}})));
end SimpleObserverStateSpaceControl;
