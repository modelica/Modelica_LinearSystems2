within Modelica_LinearSystems2.Math.Vectors;
function printStringVectorInHtml
  "Print a string vector in html format on file (without html/body heading)"
  import Modelica.Utilities.Strings;
  import Modelica.Utilities.Streams.print;

  input String s[:] "String vector";
  input String name="" "Vector name used for printing";
  input String fileName="log.html";
  input Boolean printIndices=true
    "=true, if row indices shall be printed, otherwise they are not printed";
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
        print("    <tr style=\"background-color:white\">", fileName);
        if printIndices then
           print("      <td align=\"right\" style=\"background-color:rgb(230, 230, 230);\">" + String(i) + "</td>", fileName);
        end if;
        print("      <td align=\"left\">" + s[i] + "</td>\n    </tr>", fileName);
     end for;

     // Print row closing tags
     print("    </table>", fileName);
     if name <> "" then
        print("  </td>\n</tr>\n</table>", fileName);
     end if;
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.<b>printStringVectorInHtml</b>(s, name, fileName, printIndices);
</pre></blockquote>

<h4>Description</h4>
<p>
Print vector of strings <code>s</code> called <code>name</code> into
a file <code>fileName</code> in HTML format.
Optionally, row indices are printed as well.
</p>

<h4>Example</h4>
<blockquote><pre>
s = {\"w\", \"3\", \"alpha\"};
<b>printStringVectorInHtml</b>(s, \"myVec\", \"log.html\", true);
</pre></blockquote>
<p>
saves following HTML code in the file <em>log.html</em>:
</p>
<blockquote><pre>
&lt;table border=\"0\"&gt;
&lt;tr&gt;
  &lt;td valign=\"middle\"&gt;&amp;nbsp;&amp;nbsp;&amp;nbsp;myVec&lt;/td&gt;
  &lt;td valign=\"middle\"&gt;=&lt;/td&gt;
  &lt;td&gt;
    &lt;table style=\"background-color:rgb(100, 100, 100);\" cellpadding=\"3\" border=\"0\" cellspacing=\"1\"&gt;
    &lt;tr style=\"background-color:white\"&gt;
      &lt;td align=\"right\" style=\"background-color:rgb(230, 230, 230);\"&gt;1&lt;/td&gt;
      &lt;td align=\"left\"&gt;w&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr style=\"background-color:white\"&gt;
      &lt;td align=\"right\" style=\"background-color:rgb(230, 230, 230);\"&gt;2&lt;/td&gt;
      &lt;td align=\"left\"&gt;3&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr style=\"background-color:white\"&gt;
      &lt;td align=\"right\" style=\"background-color:rgb(230, 230, 230);\"&gt;3&lt;/td&gt;
      &lt;td align=\"left\"&gt;alpha&lt;/td&gt;
    &lt;/tr&gt;
    &lt;/table&gt;
  &lt;/td&gt;
&lt;/tr&gt;
&lt;/table&gt;
</pre></blockquote>
</html>", revisions="<html>
</html>"));
end printStringVectorInHtml;
