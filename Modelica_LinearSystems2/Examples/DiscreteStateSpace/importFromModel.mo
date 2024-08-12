within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function importFromModel
  "Generate a linear discrete state space system from a (nonlinear) Modelica model"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.DiscreteStateSpace;

  input String modelName="Modelica_LinearSystems2.Utilities.Plants.DoublePendulum";
  input Modelica.Units.SI.Time Ts=0.01;
  input Modelica.Units.SI.Time T_linearize=0
    "Point in time of simulation to linearize the model";

  output DiscreteStateSpace dss=DiscreteStateSpace.Import.fromModel(
    modelName=modelName,
    Ts=Ts,
    T_linearize=T_linearize);
algorithm

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This example demonstrates how to generate a&nbsp;linear discrete state space system from
a&nbsp;(nonlinear) Modelica model
<a href=\"modelica:/Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\">DoublePendulum</a>
by its linearization.
</p>
</html>"));
end importFromModel;
