within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO1Log

algorithm
   Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
        "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO1",
        modelParam={Modelica_LinearSystems2.Records.ParameterVariation(
                     Name="k", grid=Modelica_LinearSystems2.Types.Grid.Logarithmic, nPoints=100, Min=0, Max=1000)},
        diagram=Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram(
                    linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
                    lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None));
end rootLocusOfControlledSISO1Log;
