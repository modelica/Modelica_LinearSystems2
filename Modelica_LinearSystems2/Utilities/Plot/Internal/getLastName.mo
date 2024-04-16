within Modelica_LinearSystems2.Utilities.Plot.Internal;
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
     (<a href=\"https://www.dlr.de/sr/en/\">DLR Institute of System Dynamics and Control</a>)<br>  
     The research leading to these results has received funding from the European Union’s Seventh
     Framework Programme (FP7/2007-2016) for the Clean Sky Joint Technology Initiative under
     grant agreement no. CSJU-GAM-SGO-2008-001.</td></tr>
</table>
</html>"));
end getLastName;
