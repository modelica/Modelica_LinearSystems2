within Modelica_LinearSystems2.Math.Vectors;
function householderVector
  "Calculate a normalized householder vector for the reflexion of vector a onto vector b "

  import Modelica_LinearSystems2.Math.Vectors.length;

  input Real a[:];
  input Real b[size(a, 1)];
  output Real u[size(a, 1)];
protected
  Real length_a=length(a);
  Real alpha=if length(a + length_a*b) > 0 then length_a else -length_a;

algorithm
  assert(length(b) > 0,
    "vector b in function housholderVector is zero vector, but at least one element should be different from zero ");
  assert(length(a) > 0,
    "vector a in function housholderVector is zero vector, but at least one element should be different from zero");
  u := (a + alpha*b)/length(a + alpha*b);

  annotation (Documentation(info="<HTML>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
Vectors.<b>householderVector</b>(a, b);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The function call \"<code>householderVector(a, b)</code>\" returns vector
<b>u</b>, which is the normalized Householder vector for a Householder
reflexion with matrix Q
<p>
Q = I - 2*u*u'
<p>
with
<p>
Q*a = c*b
</HTML>", revisions="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.Utilities<b>householderVector</b>(a,b);
</pre></blockquote>
<h4>Description</h4>
<p>
The function call \"<code>householderVector(a, b)</code>\" returns vector
<b>u</b>, which is the normalized Householder vector for a Householder
reflexion with matrix <b>Q</b>
</p>
<blockquote>
<p>
<b>Q</b> = <b>I</b> - 2*<b>u</b>*<b>u</b>',
</p>
</blockquote>
i.e., vector <b>a</b> is mapped to
<blockquote>
<p>
<b>a</b> -> <b>Q</b>*<b>a</b>=c*<b>b</b>
</p>
</blockquote>
with scalar c. <b>Q</b>*<b>a</b> is the reflection of <b>a</b> about the hyperplane orthogonal to <b>u</b>.
<b>Q</b> is an orthogonal matrix, i.e.
<blockquote>
<p>
    <b>Q</b> = inv(<b>Q</b>) = <b>Q</b>'
</p>
</blockquote>
<h4>Example</h4>
<blockquote><pre>
  a = {2, -4, -2, -1};
  b = {1, 0, 0, 0};

  u=<b>householderVector</b>(a,b);    // {0.837, -0.478, -0.239, -0.119}
                                      // Computation (I - 2*matrix(u)*transpose(matrix(u))) results in         
                                      // {-5, 0, 0, 0} = -5*b        
</pre></blockquote>
<h4>See also</h4>
<a href=\"modelica://Modelica_LinearSystems2.Math.Vectors.householderReflexion\">Vectors.householderReflexion</a>
</html>"));
end householderVector;
