within Modelica_LinearSystems2.Math.Vectors;
function householderVector
  "calculate a normalized householder vector for the reflexion of vector a onto vector b "

  import Modelica_LinearSystems2.Math.Vectors.length;

  input Real a[:];
  input Real b[size(a, 1)];
  output Real u[size(a, 1)];
protected
  Real la=length(a);
  Real alpha=if length(a + la*b) > 0 then la else -la;

algorithm
  assert(length(b) > 0,
    "vector b in function housholderVector is zero vector, but al least one element should be different to zero ");
  assert(length(a) > 0,
    "vector a in function housholderVector is zero vector, but al least one element should be different to zero");
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
</HTML>"));
end householderVector;
