within Modelica_LinearSystems2.Utilities.Plot.Examples;
function showMatrixDiagrams
  "Demonstrate the layout of diagrams in matrix layout"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Utilities.Plot.Records;
  import Modelica_LinearSystems2.Utilities.Plot.Types;

protected
  Real x[2]={0,1};
  Real y1[2]={0,0.5};
  Real y2[2]={0.2,0.7};
  Real y3[2]={0.4,0.9};
  Records.Curve curves[3]={
    Records.Curve(x=x, y=y1, legend="curve 1"),
    Records.Curve(x=x, y=y2, legend="curve 2"),
    Records.Curve(x=x, y=y3, legend="curve 3")};
algorithm
  Modelica_LinearSystems2.Utilities.Plot.diagramMatrix(
    fill(Records.Diagram(curve=curves), 5, 3), Records.Device(diagramWidth=80));

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots a following set of diagrams using a matrix layout (with default input arguments):
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showMatrixDiagrams.png\">
</blockquote>
</html>"));
end showMatrixDiagrams;
