within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfDriveSolidLine
  "Plot the root locus of a drive with varying load"
  extends Modelica.Icons.Function;

algorithm
  Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
    "Modelica.Mechanics.Rotational.Examples.First",
    modelParam={
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="Jload",
        grid=Modelica_LinearSystems2.Utilities.Types.Grid.Equidistant,
        Min=1,
        Max=20,
        nPoints=30),
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="Jmotor",
        Value=0.1)},
    diagram=Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram(
      linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
      lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None));

  annotation (
    Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica.Mechanics.Rotational.Examples.First\">Rotational.Examples.First</a>
over the load inertia <b>Jload</b>:
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfDrive.png\">
</blockquote>
</html>"));
end rootLocusOfDriveSolidLine;
