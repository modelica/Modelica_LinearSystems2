within Modelica_LinearSystems2.Controller.Examples;
model SimpleControlledDrive
  "Simple P-PI cascade controller to control a flexible drive"
  extends Modelica.Icons.Example;
  parameter Real kp = 10 "Gain of P position controller";
  parameter Real kv = 9 "Gain of PI speed controller";
  parameter Modelica.SIunits.Time Tv = 0.05
    "Time constant of PI speed controller";
  Modelica.Mechanics.Rotational.Components.Inertia motorInertia(J=0.1)
    annotation (extent=[-20,-80; 0,-60], Placement(transformation(extent={{-20,
            -80},{0,-60}})));
  Modelica.Mechanics.Rotational.Components.Inertia loadInertia(J=0.3)
    annotation (extent=[60,-80; 80,-60]);
  Modelica.Mechanics.Rotational.Components.SpringDamper spring(c=1e5, d=100)
    annotation (extent=[20,-80; 42,-60]);

  Modelica.Mechanics.Rotational.Sources.Torque torque
    annotation (extent=[-60,-80; -40,-60]);
  inner SampleClock sampleClock(sampleTime=0.005, blockType=
        Modelica_LinearSystems2.Controller.Types.BlockType.Discrete)
    annotation (extent=[-83,55; -63,75], Placement(transformation(extent={{-80,
            60},{-60,80}})));
  Modelica.Blocks.Sources.Ramp ramp(duration=2)
                                    annotation (extent=[-98,20; -78,40]);
  Filter filter(
    f_cut=5,
    analogFilter=Modelica_LinearSystems2.Types.AnalogFilter.Bessel,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Discrete)
             annotation (extent=[-72,20; -52,40], Placement(transformation(
          extent={{-71,20},{-51,40}})));
  Sampler sampler3(sampleFactor=2) annotation (extent=[-44,20; -24,40],
      Placement(transformation(extent={{-44,20},{-24,40}})));
  Modelica.Blocks.Math.Feedback feedback annotation (extent=[-20,20; 0,40]);
  Modelica.Blocks.Math.Gain gain(k=kp) annotation (extent=[8,20; 28,40],
      Placement(transformation(extent={{10,20},{30,40}})));
  Modelica.Mechanics.Rotational.Sensors.AngleSensor angle
    annotation (extent=[-20,-46; 0,-26], rotation=90);
  Modelica.Blocks.Math.Feedback feedback2 annotation (extent=[36,20; 56,40]);
  PI PI1(k=kv, T=Tv,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous)
                     annotation (extent=[66,20; 86,40], Placement(
        transformation(extent={{64,20},{84,40}})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speed
    annotation (extent=[36,-44; 56,-24], rotation=90);
  Sampler sampler1(sampleFactor=2)
    annotation (extent=[-20,-12; 0,8], rotation=90);
  Sampler sampler2(sampleFactor=2)
                   annotation (extent=[36,-10; 56,10], rotation=90);
equation
  connect(motorInertia.flange_b, spring.flange_a)
    annotation (points=[0,-70; 20,-70], style(color=0, rgbcolor={0,0,0}));
  connect(spring.flange_b, loadInertia.flange_a)
    annotation (points=[42,-70; 60,-70], style(color=0, rgbcolor={0,0,0}));
  connect(sampler1.y, feedback.u2) annotation (points=[-10,9; -10,22],
      style(color=74, rgbcolor={0,0,127}));
  connect(sampler2.y, feedback2.u2) annotation (points=[46,11; 46,22],
      style(color=74, rgbcolor={0,0,127}));
  connect(speed.w, sampler2.u) annotation (points=[46,-23; 46,-12], style(
        color=74, rgbcolor={0,0,127}));
  connect(angle.phi, sampler1.u) annotation (points=[-10,-25; -10,-14],
      style(color=74, rgbcolor={0,0,127}));
  connect(torque.flange, motorInertia.flange_a) annotation (Line(
      points={{-40,-70},{-20,-70}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(motorInertia.flange_a, angle.flange) annotation (Line(
      points={{-20,-70},{-26,-70},{-26,-52},{-10,-52},{-10,-46}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(motorInertia.flange_b, speed.flange) annotation (Line(
      points={{0,-70},{16,-70},{16,-50},{46,-50},{46,-44}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(filter.u, ramp.y) annotation (Line(
      points={{-73,30},{-77,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sampler3.y, feedback.u1) annotation (Line(
      points={{-23,30},{-18,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback.y, gain.u) annotation (Line(
      points={{-1,30},{8,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, feedback2.u1) annotation (Line(
      points={{31,30},{38,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback2.y, PI1.u) annotation (Line(
      points={{55,30},{62,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(torque.tau, PI1.y) annotation (Line(
      points={{-62,-70},{-82,-70},{-82,-84},{92,-84},{92,30},{85,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(filter.y, sampler3.u) annotation (Line(
      points={{-50,30},{-46,30}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(
      coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
      Rectangle(extent=[-50,60; 96,-18], style(color=1, rgbcolor={255,0,0})),
      Text(
        extent=[0,68; 60,62],
        style(color=1, rgbcolor={255,0,0}),
        string="controller"),
      Text(
        extent=[-111,12; -48,6],
        style(color=1, rgbcolor={255,0,0}),
        string="reference"),
      Text(
        extent=[-12,-91; 53,-97],
        style(color=1, rgbcolor={255,0,0}),
        string="plant (flexible drive)"),
      graphics),
    experiment(StopTime=3),
    experimentSetupOutput,
    Commands(file="Extras/Scripts/SimpleControlledDriver_Plot1.mos"
        "Plot most important variables"),
    Coordsys(grid=[1,1], scale=0),
    Documentation(info="<html>
<p>
This example demonstrates the control of a simple model
of a flexible drive system with a continuous or discrete
P-PI cascade controller. Simulate for 3 s and plot
</p>
<pre>
  ramp.y          (reference angle of loadInertia)
  loadInertia.phi (angle of loadInertia)
  loadInertia.w   (speed of loadInertia)
  torque.tau      (motor torque)         
</pre>
<p>
The standard setting in component sampleClock models a continuous controller.
This means that all 3 samplers are just dummy components containing the
equation \"y=u\" and that the PI component in the controller is a continuous
PI controller.
</p>
<p>
Change sampleClock.blockType to \"Discrete\" block. By this global setting,
the 3 sampler blocks and the PI speed controller are transformed into
a discrete representation. The base sample time is defined in
component sampleClock (= 0.02 s). Every discrete component samples
its input and output. The sampling time of every component is a multiple
of the base sample time (defined via parameter sampleFactor). 
Here, the sampler and the PI speed controller are sampled with the
base sample frequency. The sample time of the 2 samplers and
of the P position controller is a factor of 5 slower.
</p>
<p>
When comparing the simulations of the continuous and the 
(more realistic) discrete representation, it turns out that the
discrete control systems works a bit worse. This can be improved
by reducing the sample time in sampleClock.
</p>
<p>
The Controller library has several blocks to model this system
even more realistically, e.g, by component AD converter to model
the quantization errors of the analog measurement signals,
component DA converter to model the quantization errors and computing
time to determine the analog actuator (torque) signal, and
component Noise to add uniformly distributed noise to
the measurement signals.
</p>
<p>
In the following figure simulation results of the discrete and
of the continuous controller are shown:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/SimpleControlledDrive_Plot1.png\">
</p>
</html>"));
end SimpleControlledDrive;
