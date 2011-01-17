within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model disturbedSine
  import Modelica_LinearSystems2;

  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

  Real x[3];
  parameter Real a=1.0;
  Modelica_LinearSystems2.Controller.Noise w[2](y_min={-0.2,-0.1}, y_max={0.2,0.1})
    annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
  Modelica_LinearSystems2.Controller.Noise r(y_min=-0.3, y_max=0.3)
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(sampleTime=0.01,
      blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Discrete)
    annotation (Placement(transformation(extent={{60,68},{80,88}})));
equation
  der(x) = [0,1,0;0,0,0;0,0,0]*x + [0,0;1,0;0,1]*w.y;
  y = (a+x[3])*sin(x[1])+r.y;

end disturbedSine;
