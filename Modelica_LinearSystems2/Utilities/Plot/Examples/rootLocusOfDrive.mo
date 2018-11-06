within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfDrive "Plot the root locus of a drive with varying load"
  extends Modelica.Icons.Function;

algorithm
  Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
    "Modelica.Mechanics.Rotational.Examples.First",
    modelParam={
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="Jload",
        grid=Modelica_LinearSystems2.Utilities.Types.Grid.Equidistant,
        nPoints=30,
        Min=1,
        Max=20),
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="Jmotor",
        Value=0.1)});

  annotation (
    Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica.Mechanics.Rotational.Examples.First\">Rotational.Examples.First</a>
over the load inertia <b>Jload</b> (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same parameter value are
marked with a red square):
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/RootLocusOfModel.png\"/></p>
</html>"));
end rootLocusOfDrive;
