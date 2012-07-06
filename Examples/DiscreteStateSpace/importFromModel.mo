within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function importFromModel
  "This example demonstrates how to generate a linear discrete state space system from a (nonlinear) Modelica moodel"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input String modelName="Modelica_LinearSystems2.Utilities.Plants.DoublePendulum";
  input Modelica.SIunits.Time Ts=0.01;
  input Modelica.SIunits.Time T_linearize=0
    "point in time of simulation to linearize the model";

  output DiscreteStateSpace dss=DiscreteStateSpace.Import.fromModel(modelName=
    modelName, Ts=Ts,T_linearize=T_linearize);
algorithm

annotation (__Dymola_interactive=true);
end importFromModel;
