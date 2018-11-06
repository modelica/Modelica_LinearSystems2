within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO1Log "Root locus of a SISO system over controller gain k with logarithmic gridding"
  extends Modelica.Icons.Function;

algorithm
  Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
    "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO1",
    modelParam={
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="k",
        grid=Modelica_LinearSystems2.Utilities.Types.Grid.Logarithmic,
        nPoints=100,
        Min=0,
        Max=1000)},
    diagram=Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram(
      linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
      lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None));

  annotation (
    Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO1\">Plot.Examples.Utilities.ControlledSISO1</a>
over the controller gain <b>k</b> (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same parameter value are
marked with a red square):
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfControlledSISO1Log.png\"/></p>
</html>"));
end rootLocusOfControlledSISO1Log;
