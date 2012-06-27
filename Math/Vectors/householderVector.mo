within Modelica_LinearSystems2.Math.Vectors;
function householderVector
  "Calculate a normalized householder vector to reflect vector a onto vector b"

  import Modelica.Math.Vectors.norm;

  input Real a[:] "Real vector to be reflected";
  input Real b[size(a, 1)] "Real vector b vector a is mapped onto";
  output Real u[size(a, 1)] "Housholder vector to map a onto b";
protected
  Real norm_a=norm(a,2);
  Real norm_b=norm(b,2);
  Real alpha;

algorithm
  assert(norm_b > 0, "Vector b in function housholderVector is zero vector, but at least one element should be different from zero");
  assert(norm_a > 0, "Vector a in function housholderVector is zero vector, but at least one element should be different from zero");
  alpha := if norm(a + norm_a/norm_b*b,2) > norm(a - norm_a/norm_b*b,2) then norm_a/norm_b else -norm_a/norm_b;
  u := (a + alpha*b)/length(a + alpha*b);

  annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.Utilities.<b>householderVector</b>(a,b);
</pre></blockquote>
<h4>Description</h4>
<p>
The function call \"<code>householderVector(a, b)</code>\" returns the normalized Householder vector
<b>u</b> for Householder reflection of input vector <b>a</b> onto vector <b>b</b>, i.e. Householder vector <b>u</b> is the normal
vector of the reflection plane. Algebraically, the reflection is performed by transformation matrix <b>Q</b>
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
with scalar c, |c| = ||<b>a</b>|| / ||<b>b</b>||. <b>Q</b>*<b>a</b> is the reflection of <b>a</b> about the hyperplane orthogonal to <b>u</b>.
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

  u = <b>householderVector</b>(a,b);    // {0.837, -0.478, -0.239, -0.119}
                               // Computation (identity(4) - 2*matrix(u)*transpose(matrix(u)))*a results in
                               // {-5, 0, 0, 0} = -5*b


</HTML>"));
end householderVector;
