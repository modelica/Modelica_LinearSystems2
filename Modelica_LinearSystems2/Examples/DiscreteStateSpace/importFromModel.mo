within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function importFromModel
  "Example generating a discrete state space system from a model"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.DiscreteStateSpace;

  input String modelName="Modelica_LinearSystems2.Utilities.Plants.DoublePendulum";
  input Modelica.SIunits.Time Ts=0.01;
  input Modelica.SIunits.Time T_linearize=0
    "Point in time of simulation to linearize the model";

  output DiscreteStateSpace dss=DiscreteStateSpace.Import.fromModel(
    modelName=modelName,
    Ts=Ts,
    T_linearize=T_linearize);
algorithm

  annotation (
    Documentation(info="<html>
<p>
This example shows how to generate a linear discrete state space system
by linearization of a (nonlinear) Modelica model.
</p>
</html>
"));
end importFromModel;
