within Modelica_LinearSystems2.Utilities.Plot.Examples;
function rootLocusOfPIDDrive
 "Root locus of a PID controlled drive system over controller time constant Ti with logarithmic gridding"
  extends Modelica.Icons.Function;

algorithm
  Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel(
    "Modelica.Blocks.Examples.PID_Controller",
    { Modelica_LinearSystems2.Records.ParameterVariation(
        Name="PI.k",
        Value=100,
        Min=0,
        Max=1e+100),
      Modelica_LinearSystems2.Records.ParameterVariation(
        Name="PI.Ti",
        grid=Modelica_LinearSystems2.Utilities.Types.Grid.Logarithmic,
        nPoints=100,
        Value=0.1,
        Min=1e-3,
        Max=10)});

  annotation (
    Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica.Blocks.Examples.PID_Controller\">Modelica.Blocks.Examples.PID_Controller</a>
over the time constant <b>Ti</b> of the PID controller with a logaritmic gridding in Ti and a fixed
value k=100 of the PID gain k (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same parameter value are
marked with a red square):
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfPIDDrive.png\"/></p>
</html>"));
end rootLocusOfPIDDrive;
