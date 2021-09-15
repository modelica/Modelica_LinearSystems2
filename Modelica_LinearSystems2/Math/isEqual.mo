within Modelica_LinearSystems2.Math;
function isEqual
  "Determine whether two Real numbers are numerically identical (up to machine precision)"
  extends Modelica.Icons.Function;
  input Real u1 "First scalar";
  input Real u2 "Second scalar";
  input Real eps(min=0) = 0.0
    "The two scalars are identical if abs(u1-u2) <= eps";
  output Boolean result "True, if abs(u1-u2) <= eps";
algorithm
  result :=abs(u1 - u2) <= eps;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = <strong>isEqual</strong>(r1, r2);
result = <strong>isEqual</strong>(r1, r2, eps=0.0);
</pre></blockquote>

<h4>Description</h4>

<p>
The function call &quot;<code>isEqual(r1, r2)</code>&quot; returns <strong>true</strong>,
if the two Real numbers r1 and r2 are the same up to a given precision eps.
(result = abs(r1-r2) &le; eps). Otherwise the function
returns <strong>false</strong>. With the optional third argument <strong>eps</strong>
the range can be defined, in which two Real numbers are treated as identical.
The default is &quot;eps = 0&quot;. Another useful value is, e.g.,
&quot;eps = 10*Modelica.Constants.eps&quot;.
</p>

<h4>Example</h4>
<blockquote><pre>
  Real r1 = 3;
  Real r2 = 3;
  Real r3 = 3.0001;
  Boolean result;
<strong>algorithm</strong>
  result := isEqual(r1,r2);          // = <strong>true</strong>
  result := isEqual(r1,r3);          // = <strong>false</strong>
  result := isEqual(r1,r3, eps=0.1); // = <strong>true</strong>
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Vectors.isEqual\">Vectors.isEqual</a>,
<a href=\"modelica://Modelica.Math.Matrices.isEqual\">Matrices.isEqual</a>,
<a href=\"modelica://Modelica.Utilities.Strings.isEqual\">Strings.isEqual</a>
</p>
</html>"));
end isEqual;
