within Modelica_LinearSystems2.Math.Vectors;
function printStringVectorInHtml
  "Print a string vector in html format on file (without html/body heading)"
  import Modelica.Utilities.Strings;
  import Modelica.Utilities.Streams.print;

  input String s[:] "String vector";
  input String name="" "Vector name used for printing";
  input String fileName="log.html";
protected
  Integer r=size(s, 1);
  Boolean siExist=false;
algorithm
  // Check whether all entries are an empty string
  for i in 1:r loop
     siExist := siExist or (s[i] <> "");
  end for;

  if not siExist then
     if name == "" then
        print("<p>&nbsp;&nbsp;&nbsp;[&nbsp;], empty vector</p>\n", fileName);
     else
        print("<p>&nbsp;&nbsp;&nbsp;" + name + "&nbsp;=&nbsp;[&nbsp;], empty vector</p>\n", fileName);
     end if;
  else
     // Print "name = "
     if name <> "" then
        print("<table border=\"0\">\n<tr>\n  <td valign=\"middle\">" +
              "&nbsp;&nbsp;&nbsp;" + name + "</td>\n  <td valign=\"middle\">=</td>\n  <td>", fileName);
     end if;

     // Print table heading
     print("    <table style=\"background-color:rgb(100, 100, 100);\" cellpadding=\"3\" border=\"0\" cellspacing=\"1\">",
           fileName);

     // Print vector elements
     for i in 1:r loop
        print("    <tr style=\"background-color:white\">\n      <td align=\"right\" style=\"background-color:rgb(230, 230, 230);\">"
                  + String(i) + "</td>", fileName);
        print("      <td align=\"left\">" + s[i] + "</td>\n    </tr>", fileName);
     end for;

     // Print row closing tags
     print("    </table>", fileName);
     if name <> "" then
        print("  </td>\n</tr>\n</table>", fileName);
     end if;
  end if;

  annotation (Documentation(info="<html>

</html>", revisions="<html>
</html>"));
end printStringVectorInHtml;
