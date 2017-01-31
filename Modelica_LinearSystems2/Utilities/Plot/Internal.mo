within Modelica_LinearSystems2.Utilities.Plot;
package Internal "Internal functions, that should not be utilized by a user"
  extends Modelica.Icons.Package;
  function convertToDymolaMarker
    "Convert from marker type of Modelica_LinearSystems2 to Dymola marker"
    import Symbol = Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol;

    input Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol symbol
      "Enumeration value of point symbol used in Modelica_LinearSystems2";
    output MarkerStyle result
      "Enumeration value of point symbol used in Dymola";
  algorithm
    result := if symbol == Symbol.None then MarkerStyle.None else if symbol ==
      Symbol.Cross then MarkerStyle.Cross else if symbol == Symbol.Circle then
      MarkerStyle.Circle else if symbol == Symbol.Square then MarkerStyle.Square
       else if symbol == Symbol.FilledSquare then MarkerStyle.FilledSquare
       else if symbol == Symbol.TriangleDown then MarkerStyle.TriangleDown
       else if symbol == Symbol.TriangleUp then MarkerStyle.TriangleUp else if
      symbol == Symbol.Diamond then MarkerStyle.Diamond else if symbol ==
      Symbol.Dot then MarkerStyle.Dot else MarkerStyle.None;

    annotation (Inline=true, Documentation(info="<html>
<p>
This function converts the point symbol enumeration as used in the Modelica_LinearSystems2 library
to the enumeration used in Dymola.
</p>
</html>"));
  end convertToDymolaMarker;

  function convertToDymolaPattern
    "Convert from pattern type of Modelica_LinearSystems2 to Dymola pattern"
    import Pattern = Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern;

    input Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern pattern
      "Enumeration value of line pattern used in Modelica_LinearSystems2";
    output LinePattern result
      "Enumeration value of line pattern used in Dymola";
  algorithm
    result := if pattern == Pattern.None then LinePattern.None else if pattern ==
      Pattern.Solid then LinePattern.Solid else if pattern == Pattern.Dash
       then LinePattern.Dash else if pattern == Pattern.Dot then LinePattern.Dot
       else if pattern == Pattern.DashDot then LinePattern.DashDot else if
      pattern == Pattern.DashDotDot then LinePattern.DashDotDot else
      LinePattern.Solid;

    annotation (Inline=true, Documentation(info="<html>
<p>
This function converts the line pattern enumeration as used in the Modelica_LinearSystems2 library
to the enumeration used in Dymola.
</p>
</html>"));
  end convertToDymolaPattern;

  function getLastName "Get last name of a Modelica path name"
     import Modelica.Utilities.Strings;
     input String path;
     output String tail "= last part of path (after the last '.'";
  protected
     Integer startIndex;
     Integer endIndex;
  algorithm
     startIndex :=Strings.findLast(path, ".");
     if startIndex == 0 or startIndex >= Strings.length(path) then
        tail := path;
     else
        tail := Strings.substring(path, startIndex+1, Strings.length(path));
     end if;
    annotation (Documentation(revisions="<html>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td valign=\"top\"> Nov. 29, 2015 </td>
    <td valign=\"top\">
     Initial version implemented by
     Martin R. Kuhn and Martin Otter 
     (<a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>)<br>  
     The research leading to these results has received funding from the European Union’s Seventh
     Framework Programme (FP7/2007-2016) for the Clean Sky Joint Technology Initiative under
     grant agreement no. CSJU-GAM-SGO-2008-001.</td></tr>
</table>
</html>"));
  end getLastName;
  annotation (Documentation(info="<html>
<p>
This package contains utility functions to implement the plot functions.
These functions should not be used by a user.
</p>
</html>"));
end Internal;
