within Modelica_LinearSystems2.Math.Vectors;
function length "Return length of a vector"

  input Real v[:] "Vector";
  output Real result "Length of vector v";

algorithm
  result := sqrt(v*v);
  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
Vectors.<b>length</b>(v);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>Vectors.<b>length</b>(v)</code>\" returns the
<b>Euclidean length</b> \"<code>sqrt(v*v)</code>\" of vector v. 
The function call is equivalent to Vectors.norm(v). The advantage of
length(v) over norm(v)\"is that function length(..) is implemented
in one statement and therefore the function is usually automatically
inlined. Further symbolic processing is therefore possible, which is
not the case with function norm(..).
</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  v = {2, -4, -2, -1};
  <b>length</b>(v);  // = 5
</pre></blockquote>
<h4><font color=\"#008000\">See also</font></h4>
<a href=\"Modelica:Modelica_LinearSystems2.Math.Vectors.norm\">Vectors.norm</a>
</html>"));
end length;
