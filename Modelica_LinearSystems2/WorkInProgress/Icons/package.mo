within Modelica_LinearSystems2.WorkInProgress;
package Icons "Package of icons used in the Linear Systems library"
  extends Modelica.Icons.IconsPackage;


annotation (Documentation(info="<html>
<p>
  Package contains Icons used in the Linear Systems library.
</p>

<h4>Graphical display of function maturity</h4>
<p>
Mainly to facilitate the development process following icons exist to
identify the level of maturity of a function:
</p>
<ol>
<li>NotWorkingYetFunction - Functions that are not working yet or a known to having severe problems</li>
<li>DeveloperFunction - Functions that are in an early development stage. They probably already provide some usefullnes</li>
<li>Release level function - This functions use either the standard Modelica function icon, or no icon at all</li>
</ol>
<p>
Additionally following icons are provided that can be combined with the former ones:
</p>
<ul>
<li>DeprecatedFunction - Deprecated functions that should not be used any more</li>
<li>FineDocumentationFunction - Functions which pretty well documented.
As soon as there are more documented than non-documented function, this Icon may be replaced by a BadDocumentationFunction-Icon.</li>
</ul>

  </html>"));
end Icons;
