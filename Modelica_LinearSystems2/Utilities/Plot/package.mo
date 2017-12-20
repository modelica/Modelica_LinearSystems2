within Modelica_LinearSystems2.Utilities;
package Plot "Package of functions for generation of 2D-plots"
  extends Modelica.Icons.Package;

  annotation (Documentation(info="<html>
<p>
This package provides functions to plot curves in two dimensions. Here is a short overview:
</p>

<p>
A figure consists of a <b>set of diagrams</b>. Different functions are provided
to either plot one diagram (Plot.diagram) to plot several diagrams under each other
(Plot.vectorDiagrams) or to plot several diagrams in matrix layout
(Plot.matrixDiagrams).
</p>

<p>
Every diagram can have a set of <b>curves</b>. Every diagram has the same
width, defined by <b>diagramWidth</b>. The height of a diagram is defined
by variable <b>heightRatio</b>
(diagram height in row j = diagram[j].heightRatio*diagramWidth).
Several curves can be displayed in one diagram.
</p>

<p>
A typical example is shown in the next figure:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showSinesInVectorDiagrams.png\"></p>
</html>"));
end Plot;
