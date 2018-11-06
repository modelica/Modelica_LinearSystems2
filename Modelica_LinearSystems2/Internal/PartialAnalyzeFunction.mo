within Modelica_LinearSystems2.Internal;
partial function PartialAnalyzeFunction
  "Partial function to linearize a model and then perform operations on the linearized model"

  input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
  input Modelica_LinearSystems2.Records.SetParameter modelParam[:]=
    fill(Modelica_LinearSystems2.Records.SetParameter(Name="",Value=0.0),0)
    "Values of model parameters used for linearization";
  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options" annotation(Dialog(enable=not linearizeAtInitial));

protected
  Modelica_LinearSystems2.StateSpace ssLin=
       Modelica_LinearSystems2.Utilities.Import.linearize2(modelName, modelParam, simulationSetup);
end PartialAnalyzeFunction;
