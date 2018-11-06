within Modelica_LinearSystems2.Utilities.Plot;
function diagramMatrix "Plot several diagrams in matrix layout"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram[:,:]
    "Properties of a set of diagrams (matrix layout)"
    annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
    Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);
protected
  Modelica_LinearSystems2.Utilities.Plot.Records.Device device2=device;

algorithm
  for i in 1:size(diagram,2) loop
    device2.xTopLeft :=device.xTopLeft + (i - 1)*device.diagramWidth;
    Modelica_LinearSystems2.Utilities.Plot.diagramVector(diagram[:, i], device2);
  end for;

  annotation (
    Documentation(info="<html>
<p>
This function plots a set of 2-dimensional curves in a set of diagrams
using a matrix layout. For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>

<h4>Example</h4>
<p>
See an <a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.showMatrixDiagrams\">example</a>
for possible usage of this function:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showMatrixDiagrams.png\"></p>
</html>"));
end diagramMatrix;
