within Modelica_LinearSystems2.Utilities;
package Plot "Functions for generation of 2D-plots"
    extends Modelica.Icons.Package;






  annotation (Documentation(info="<html>
<p>
This package provides functions to plot curves in two dimensions. Here is a short overview:
</p>

<p>
A figure consists of a <u>set of diagrams</u>. Different functions are provided
to either plot one diagram (Plot.diagram) to plot several diagrams under each other
(Plot.vectorDiagrams) or to plot several diagrams in matrix layout
(Plot.matrixDiagrams).
</p>

<p>
Every diagram can have a set of <u>curves</u>. Every diagram has the same
width, defined by <u>diagramWidth</u>. The height of a diagram is defined
by variable <u>heightRatio</u> 
(diagram height in row j = diagram[j].heightRatio*diagramWidth). 
Several curves can be displayed in one diagram.
</p>
 
</html>"));
end Plot;
