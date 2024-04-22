within Modelica_LinearSystems2.Controllers.BaseClasses;
partial model SimpleObserverStateSpaceControl
  "Template for a simple state feedback controller with observer and optional pre-filter"

  MatrixGain feedbackMatrix
    annotation (Placement(transformation(extent={{30,-50},{10,-30}})));
  MatrixGain preFilter(K=[0])
    annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
  Modelica.Blocks.Math.Feedback feedback[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  replaceable PartialPlantMIMO plant(
    n=feedbackMatrix.nin,
    m=feedbackMatrix.nout) constrainedby PartialPlantMIMO
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  Sampler samplerPreFilter[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{-43,-5},{-33,5}})));
  Sampler samplerFeedback[feedbackMatrix.nout]
    annotation (Placement(transformation(extent={{-5,-45},{-15,-35}})));
  Sampler samplerOut[observer.nout]
    annotation (Placement(transformation(extent={{85,-51},{75,-41}})));
  inner SampleClock sampleClock
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
  Modelica_LinearSystems2.Controllers.Observer observer
    annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
equation
  connect(feedback.y, plant.u) annotation (Line(
      points={{-11,0},{78,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.u, feedbackMatrix.y) annotation (Line(
      points={{-4,-40},{9,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerFeedback.y, feedback.u2) annotation (Line(
      points={{-15.5,-40},{-20,-40},{-20,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preFilter.y, samplerPreFilter.u) annotation (Line(
      points={{-49,0},{-44,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerPreFilter.y, feedback.u1) annotation (Line(
      points={{-32.5,0},{-28,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(samplerOut.y, observer.y) annotation (Line(
      points={{74.5,-46},{62,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(plant.ym, samplerOut.u) annotation (Line(
      points={{90,-11},{90,-46},{86,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(observer.x_estimated, feedbackMatrix.u) annotation (Line(
      points={{39,-40},{32,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(observer.u, feedback.y) annotation (Line(points={{62,-34},{70,-34},{
          70,0},{-11,0}}, color={0,0,127}));
  annotation (    Documentation(info="<html>
<p>This template represents the structure of a simple state feedback controller with observer and optional pre-filter. In the application it must be extended with an input signal.</p>
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})));
end SimpleObserverStateSpaceControl;
