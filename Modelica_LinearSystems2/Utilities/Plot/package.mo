within Modelica_LinearSystems2.Utilities;
package Plot "Package of functions for generation of 2D-plots"
  extends Modelica.Icons.Package;

  annotation (Documentation(info="<html>
<p>
This package provides functions to plot curves in two dimensions. Here is a short overview:
</p>

<p>
A figure consists of a <strong>set of diagrams</strong>. Different functions are provided
to either plot one diagram (<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.diagram\">Plot.diagram</a>) to plot several diagrams under each other
(<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.diagramVector\">Plot.diagramVector</a>) or to plot several diagrams in matrix layout
(<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.diagramMatrix\">Plot.diagramMatrix</a>).
</p>

<p>
Every diagram can have a set of <strong>curves</strong>. Every diagram has the same
width, defined by <strong>diagramWidth</strong>. The height of a diagram is defined
by variable <strong>heightRatio</strong>
(diagram height in row j = diagram[j].heightRatio*diagramWidth).
Several curves can be displayed in one diagram.
</p>

<p>
A typical example is shown in the next figure:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showSinesInVectorDiagram.png\"></p>
</html>"));
end Plot;
