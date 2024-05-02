within Modelica_LinearSystems2.Math.Matrices;
function printMatrixInHtml
  "Print a matrix in html format on file (without html/body heading)"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;

  input Real M[:,:] "Real matrix";
  input String name="" "Matrix name used for printing";
  input String fileName="log.html"
    "Name of file to be printed in (incl. file extension)";
  input String format=".6g" "Format of numbers (e.g. \"20.8e\")";
  input Boolean printIndices=true
    "=true, if row and column indices shall be printed, otherwise they are not printed";
protected
  Integer r=size(M, 1);
  Integer c=size(M, 2);

algorithm
  if r == 0 or c == 0 then
     if name == "" then
        print("<p>&nbsp;&nbsp;&nbsp;[&nbsp;], empty matrix: " + String(r)+" by "+String(c) + "</p>\n", fileName);
     else
        print("<p>&nbsp;&nbsp;&nbsp;" + name + "&nbsp;=&nbsp;[&nbsp;], empty matrix: " + String(r)+" by "+String(c) + "</p>\n", fileName);
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
     if printIndices then
        print(  "\n    <tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
              + "\n      <td> </td>", fileName);

        for j in 1:c loop
           print("      <td align=\"center\">" + String(j) + "</td>", fileName);
        end for;
        print("    </tr>", fileName);
     end if;

     // Print matrix elements
     for i in 1:r loop
        print("    <tr style=\"background-color:white\">\n", fileName);
        if printIndices then
           print("      <td align=\"right\" style=\"background-color:rgb(230, 230, 230);\">"
                 + String(i) + "</td>", fileName);
        end if;

        for j in 1:c loop
           print("      <td align=\"right\">" + String(M[i, j],format=format) + "</td>", fileName);
        end for;

        print("    </tr>", fileName);
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
end printMatrixInHtml;
