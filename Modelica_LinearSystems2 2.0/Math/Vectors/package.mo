within Modelica_LinearSystems2.Math;
package Vectors "Functions operating on vectors"
extends Modelica.Icons.Library;


annotation (Documentation(info="<HTML>
<h4><font color=\"#008000\">Library content</font></h4>
<p>
This library provides functions operating on vectors:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><th><i>Function</i></th>
      <th><i>Description</i></th>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.isEqual\">isEqual</a>(v1, v2)</td>
      <td>Determines whether two vectors have the same size and elements</td>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.norm\">norm</a>(v,p)</td>
      <td>p-norm of vector v</td>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.length\">length</a>(v)</td>
      <td>Length of vector v (= norm(v,2), but inlined and therefore usable in
          symbolic manipulations) </td>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.normalize\">normalize</a>(v)</td>
      <td>Return normalized vector such that length = 1 and prevent
          zero-division for zero vector</td>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.reverse\">reverse</a>(v)</td>
      <td>Reverse vector elements</td>
  </tr>
  <tr><td><a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.sort\">sort</a>(v)</td>
      <td>Sort elements of vector in ascending or descending order</td>
  </tr>
</table>
<h4><font color=\"#008000\">See also</font></h4>
<a href=\"Modelica:Modelica.Math.Matrices\">Matrices</a>
</HTML>"));

end Vectors;
