within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function importFromModel
  "This example demonstrates how to generate a linear discrete state space system from a (nonlinear) Modelica moodel"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input String modelName="Modelica_LinearSystems2.Controller.Examples.Components.DoublePendulum";
  input Modelica.SIunits.Time Ts=0.01;

  output DiscreteStateSpace dss=DiscreteStateSpace.Import.fromModel(modelName=
    modelName, Ts=Ts);
algorithm
annotation (interactive=true);
end importFromModel;
