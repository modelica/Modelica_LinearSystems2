within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO1Log

algorithm
   Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
        "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO1",
        modelParam={Modelica_LinearSystems2.Records.ParameterVariation(
                     Name="k", nVar=1000, logVar=true, Value=1, Min=0, Max=1000)},
        diagram=Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram(
                    linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
                    lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None));
end rootLocusOfControlledSISO1Log;
