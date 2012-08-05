within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO2Log

algorithm
   Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
        "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO2",
        modelParam={Modelica_LinearSystems2.Records.ParameterVariation(
                      Name="k", nVar=100, logVar=true, Min=0.0001, Max=10)});
end rootLocusOfControlledSISO2Log;
