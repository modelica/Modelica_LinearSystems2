within Modelica_LinearSystems2.WorkInProgress.Tests.temp;
block blockFunction

  Modelica.Blocks.Interfaces.RealInput u1[2]
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealInput u2[2]
    annotation (Placement(transformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealOutput y1[size(u1,1)]
    annotation (Placement(transformation(extent={{100,30},{120,50}})));
  Modelica.Blocks.Interfaces.RealOutput y2[size(u1,1)]
    annotation (Placement(transformation(extent={{100,-70},{120,-50}})));

equation
  (y1,y2)=Modelica_LinearSystems2.WorkInProgress.Tests.temp.endFunction(function Modelica_LinearSystems2.WorkInProgress.Tests.temp.instFunction(),u1,u2);
  annotation (Diagram(graphics));
end blockFunction;
