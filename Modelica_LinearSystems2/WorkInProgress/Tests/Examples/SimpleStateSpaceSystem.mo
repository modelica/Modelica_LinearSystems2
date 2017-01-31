within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model SimpleStateSpaceSystem
  Modelica.Blocks.Continuous.StateSpace stateSpace(
    A=[2, -3; 4, 8],
    B=[3; 9],
    C=[-3, 5],
    D=[3]) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  Modelica.Blocks.Interfaces.RealInput u[size(stateSpace.u, 1)]
    "Connector of Real input signals"
    annotation (Placement(transformation(extent={{-100,30},{-60,70}})));
  Modelica.Blocks.Interfaces.RealOutput y[size(stateSpace.y, 1)]
    "Connector of Real output signals"
    annotation (Placement(transformation(extent={{8,40},{28,60}})));
equation
  connect(stateSpace.u, u) annotation (Line(
      points={{-42,50},{-80,50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(stateSpace.y, y) annotation (Line(
      points={{-19,50},{18,50}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end SimpleStateSpaceSystem;
