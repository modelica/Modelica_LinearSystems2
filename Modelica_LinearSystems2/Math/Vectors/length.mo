within Modelica_LinearSystems2.Math.Vectors;
function length
  "Return length of a vector (inlined and therefore usable in symbolic manipulations)"

  input Real v[:] "Vector";
  output Real result "Length of vector v";

algorithm
  result := sqrt(v*v);
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Vectors.length instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.<strong>length</strong>(v);
</pre></blockquote>

<h4>Description</h4>
<p>
The function call \"<code>Vectors.<strong>length</strong>(v)</code>\" returns the
<strong>Euclidean length</strong> \"<code>sqrt(v*v)</code>\" of vector v.
The function call is equivalent to <a href=\"Modelica://Modelica.Math.Vectors.norm\">Modelica.Math.Vectors.norm(v)</a>. The advantage of
length(v) over norm(v) is that function length(..) is implemented
in one statement and therefore the function is usually automatically
inlined. Further symbolic processing is therefore possible, which is
not the case with function norm(..).
</p>

<h4>Example</h4>
<blockquote><pre>
v = {2, -4, -2, -1};
<strong>length</strong>(v);  // = 5
</pre></blockquote>
</html>"));
end length;
