within Modelica_LinearSystems2.Utilities.Plot.Types;
type LinePattern = enumeration(
    None,
    Solid,
    Dash,
    Dot,
    DashDot,
    DashDotDot) "Choices for line pattern" annotation (Documentation(info="<html>
<p>
Enumeration to define the line pattern, that is the line style how the
defined points are connected together by a polyline. Possible values:
</p>

<table border=1 cellspacing=0 cellpadding=2>
<tr><th><strong>Types.LinePattern.</strong></th><th><strong>Meaning</strong></th></tr>
<tr><td valign=\"top\">None</td>
    <td valign=\"top\">Points are not connected</td></tr>

<tr><td valign=\"top\">Solid</td>
    <td valign=\"top\">Points are connected with a solid line</td></tr>

<tr><td valign=\"top\">Dash</td>
    <td valign=\"top\">Points are connected with a dash line</td></tr>

<tr><td valign=\"top\">Dot</td>
    <td valign=\"top\">Points are connected with a dotted line</td></tr>

<tr><td valign=\"top\">DashDot</td>
    <td valign=\"top\">Points are connected with a dash-dotted line</td></tr>

<tr><td valign=\"top\">DashDotDot</td>
    <td valign=\"top\">Points are connected with a dash-dotted-dotted line</td></tr>
</table>

<p>
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.showLinePatterns\">Example</a>:
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showLinePatterns.png\"></p>
</html>"));
