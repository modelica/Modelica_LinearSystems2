within Modelica_LinearSystems2.Math.Matrices.Internal;
function hessenberg2
  "Compute an upper Hessenberg matrix by repeatedly applicated householder similarity transformation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real H[:,:];
  input String s="u";
  output Real Ht[size(H, 1),size(H, 2)];

protected
  Integer q=size(H, 1);
  Real u[:] "householder vector";
  Integer ll;

algorithm
  assert(s == "u" or s == "l",
    "parameter s should be equal to 'u' or 'l' to indicate upper or lower Hessenberg form");

  Ht := H;

  for ll in 1:q - 2 loop
    u := if s == "u" then cat(
      1,
      zeros(ll),
      cat(1, Vectors.householderVector(vector(Ht[ll + 1:q, ll]), cat(
        1,
        {1},
        zeros(q - ll - 1))))) else cat(
      1,
      cat(1, Vectors.householderVector(vector(Ht[1:q - ll, q - ll + 1]), cat(
        1,
        zeros(q - ll - 1),
        {1}))),
      zeros(ll));

    Ht := if s == "u" then
      Matrices.Internal.hohoTrafoUpperHess(
      Ht,
      u,
      ll) else
      Matrices.Internal.hohoTrafoLowerHess(
      Ht,
      u,
      ll);

  end for;

  annotation (Documentation(info="<html>
<p>
This function computes the Hessenberg matrix of matrix <strong>A</strong> by repetitive application of Householder similarity transformation
</p>
<pre>
    <strong>A</strong>i+1 = (<strong>I</strong>-2*<strong>u</strong>_i*<strong>u</strong>_i')*<strong>A</strong>i*(<strong>I</strong>-2*<strong>u</strong>_i*<strong>u</strong>_i')
</pre>
<p>
with Householder vector <strong>u</strong>_i.
</p>
<p>
The elementary transformations can be subsumed under
</p>
<pre>
   <strong>A</strong> -> <strong>Q</strong>*<strong>A</strong>*<strong>Q</strong>
</pre>
<p>
and
<strong>Q</strong>*<strong>A</strong>*<strong>Q</strong> is Hessenberg matrix.
</p>
<p>
In contrast to function <strong>Modelica_LinearSystems2.Math.Matrices.hess</strong>, function <strong>Modelica_LinearSystems2.Math.Matrices.hess3</strong> does not use any LAPACK routine.
</p>
</html>"));
end hessenberg2;
