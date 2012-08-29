within Modelica_LinearSystems2.Internal;
partial function PartialAnalyzeFunction
  "Partial function to linearize a model and then perform operations on the linearized model"

  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Boolean linearizeAtInitial=true
    "= true, if linearization at inital time; otherwise simulate until t_linearize"
     annotation (choices(__Dymola_checkBox=true));
  input Modelica.SIunits.Time t_linearize= 0
    "Simulate until t_linearize and then linearize, if linearizeAtInitial=false"
                            annotation(Dialog(enable=not linearizeAtInitial));
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options if linearizeAtInitial=false" annotation(Dialog(enable=not linearizeAtInitial));
protected
  Modelica_LinearSystems2.StateSpace ssLin=
       Modelica_LinearSystems2.Utilities.Import.linearize2(modelName, linearizeAtInitial, t_linearize, simulationSetup);
annotation(__Dymola_interactive=true);
end PartialAnalyzeFunction;
