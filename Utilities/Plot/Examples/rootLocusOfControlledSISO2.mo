within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO2

algorithm
   Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
        "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO2",
        modelParam=  {Modelica_LinearSystems2.Records.ParameterVariation(Name="k", nVar=100, Min=0, Max=10)});
end rootLocusOfControlledSISO2;
