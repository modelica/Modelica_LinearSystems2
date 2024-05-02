within Modelica_LinearSystems2.Math.Vectors;
function householderVector
  "Calculate a normalized householder vector to reflect vector a onto vector b"
  extends Modelica.Icons.Function;

  import Modelica.Math.Vectors.norm;
  import Modelica.Math.Vectors.length;

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

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.<strong>householderVector</strong>(a,b);
</pre></blockquote>

<h4>Description</h4>
<p>
The function call \"<code>householderVector(a, b)</code>\" returns the normalized Householder vector
<strong>u</strong> for Householder reflection of input vector <strong>a</strong> onto vector <strong>b</strong>, i.e. Householder vector <strong>u</strong> is the normal
vector of the reflection plane. Algebraically, the reflection is performed by transformation matrix <strong>Q</strong>
</p>
<blockquote>
  <strong>Q</strong> = <strong>I</strong> - 2*<strong>u</strong>*<strong>u</strong>',
</blockquote>
<p>
i.e., vector <strong>a</strong> is mapped to
</p>
<blockquote>
  <strong>a</strong> -> <strong>Q</strong>*<strong>a</strong>=c*<strong>b</strong>
</blockquote>
<p>
with scalar c, |c| = ||<strong>a</strong>|| / ||<strong>b</strong>||. <strong>Q</strong>*<strong>a</strong> is the reflection of <strong>a</strong> about the hyperplane orthogonal to <strong>u</strong>.
<strong>Q</strong> is an orthogonal matrix, i.e.
</p>
<blockquote>
<strong>Q</strong> = inv(<strong>Q</strong>) = <strong>Q</strong>'.
</blockquote>

<h4>Example</h4>
<blockquote><pre>
a = {2, -4, -2, -1};
b = {1, 0, 0, 0};

u = <strong>householderVector</strong>(a,b);
// {0.837, -0.478, -0.239, -0.119}
// Computation (identity(4) - 2*matrix(u)*transpose(matrix(u)))*a results in
// {-5, 0, 0, 0} = -5*b
</pre></blockquote>
</html>"));
end householderVector;
