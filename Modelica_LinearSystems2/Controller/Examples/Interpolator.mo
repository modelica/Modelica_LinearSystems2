within Modelica_LinearSystems2.Controller.Examples;
model Interpolator "Demonstrate usage of Interpolator"
  extends Modelica.Icons.Example;
  import Modelica_LinearSystems2;

  Modelica.Blocks.Sources.Sine sine(
    f=2,
    offset=0.1,
    startTime=0.1)
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
  inner SampleClock sampleClock(sampleTime=0.001, blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica_LinearSystems2.Controller.Interpolator interpolator1(
    outputSampleFactor=4,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous,
    meanValueFilter=false,
    inputSampleFactor=20)
    annotation (Placement(transformation(extent={{0,40},{20,60}})));

  Modelica_LinearSystems2.Controller.Interpolator interpolator2(
    outputSampleFactor=4,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    meanValueFilter=false,
    inputSampleFactor=20)
    annotation (Placement(transformation(extent={{0,0},{20,20}})));

  Modelica_LinearSystems2.Controller.Interpolator interpolator3(
    outputSampleFactor=4,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete,
    inputSampleFactor=20)
    annotation (Placement(transformation(extent={{0,-40},{20,-20}})));
  Sampler sampler(sampleFactor=20)
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

equation
  connect(sine.y, sampler.u) annotation (Line(
      points={{-59,10},{-42,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.y, interpolator2.u)
    annotation (Line(
      points={{-19,10},{-2,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, interpolator1.u) annotation (Line(
      points={{-59,10},{-50,10},{-50,50},{-2,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler.y, interpolator3.u) annotation (Line(
      points={{-19,10},{-10,10},{-10,-30},{-2,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation ( Documentation(info="<html>
<p>
This example demonstrates the usage of the
<a href=\"modelica://Modelica_LinearSystems2.Controller.Interpolator\">Interpolator</a>
block. This block is used in multi-rate controllers to interpolate between two different
sample rates. In this example, a sampled sine-signal is interpolated in different ways:
</p>

<ol>
<li> A &quot;sine&quot; signal (= blue and red curve in the figure below) is sampled with 20*Ts
     (= green curve in the figure below) where Ts is the sample time of 1 ms.</li>

<li> This &quot;sine&quot; signal is interpolated with &quot;interpolator1&quot; which has
     a &quot;continous&quot; blockType. In this case, the interpolator just passes the input
     signal (y = u) and therefore the signal is identical to the sine signal
     (therefore the blue and the red curve are the same).</li>

<li> The sampled sine signal is interpolated  with &quot;interpolator2&quot; to arrive
     at a sample time of 4*Ts (= violet curve in the figure below).
     This gives a 20*Ts delay and a sine signal
     sampled with 4*Ts (instead of 20*Ts).</li>

<li> In order to remove frequencies introduced by the interpolation, &quot;interpolator3&quot;
     additionally filters the linearly interpolated signal with a mean value filter
     of the same length (= black curve in the figure below).
     This gives an additional delay (= a total delay of 24*Ts), but
     results in a filtered interpolation signal.</li>
</ol>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/Examples/Interpolator.png\">
</p>
</html>"),
    experiment(StopTime=0.5));
end Interpolator;
