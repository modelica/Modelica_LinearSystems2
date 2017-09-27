within Modelica_LinearSystems2.Utilities.Plot.Types;
type PointSymbol = enumeration(
    None "no symbol",
    Cross "x",
    Circle "o",
    Square,
    FilledSquare,
    TriangleDown,
    TriangleUp,
    Diamond,
    Dot) "Choices for point symbol" annotation (Documentation(info="<html>
<p>
Enumeration to define the style of the symbol to mark a point in a diagram. Possible values:
</p>

<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.PointSymbol.</b></th><th><b>Meaning</b></th></tr>
<tr><td valign=\"top\">None</td>
    <td valign=\"top\">No symbol (point is not explicitly marked)</td></tr>

<tr><td valign=\"top\">Cross</td>
    <td valign=\"top\">x</td></tr>

<tr><td valign=\"top\">Circle</td>
    <td valign=\"top\">o</td></tr>

<tr><td valign=\"top\">Square</td>
    <td valign=\"top\">Non-filled square</td></tr>

<tr><td valign=\"top\">FilledSquare</td>
    <td valign=\"top\">Filled square</td></tr>

<tr><td valign=\"top\">TriangleDown</td>
    <td valign=\"top\">Filled triangle pointing downwards</td></tr>

<tr><td valign=\"top\">TriangleUp</td>
    <td valign=\"top\">Filled triangle pointing upwards</td></tr>

<tr><td valign=\"top\">Diamond</td>
    <td valign=\"top\">Filled diamond</td></tr>

<tr><td valign=\"top\">Dot</td>
    <td valign=\"top\">.</td></tr>
</table>


<p>
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.showPointSymbols\">Example</a>:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showPointSymbols.png\"></p>
</html>"));
