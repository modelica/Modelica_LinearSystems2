within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfControlledSISO2 "Root locus of a SISO system over controller gain k with equidistant gridding"
  extends Modelica.Icons.Function;

algorithm
  Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
    "Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO2",
    modelParam={
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="k",
        grid=Modelica_LinearSystems2.Utilities.Types.Grid.Equidistant,
        nPoints=200,
        Min=0,
        Max=10)});

  annotation (
    Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.Utilities.ControlledSISO2\">Plot.Examples.Utilities.ControlledSISO2</a>
over the controller gain <b>k</b> with an equidistant gridding in k (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same parameter value are
marked with a red square):
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfControlledSISO2.png\"/></p>

<p>
Compare this plot with the logarithmic gridding in k in example
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfControlledSISO2Log\">rootLocusOfControlledSISO2Log</a>.
As can be seen, the logarithmic gridding over a controller gain yields a much better result as an equidistant gridding.
</p>
</html>"));
end rootLocusOfControlledSISO2;
