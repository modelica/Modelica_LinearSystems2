within Modelica_LinearSystems2.Math.Matrices.Internal;
function hessenberg2
  "Compute an upper Hessenberg matrix by repeatedly applicated householder similarity transformation"
  import Modelica;
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
      Modelica_LinearSystems2.Math.Matrices.Internal.hohoTrafoUpperHess(
      Ht,
      u,
      ll) else
      Modelica_LinearSystems2.Math.Matrices.Internal.hohoTrafoLowerHess(
      Ht,
      u,
      ll);

  end for;

  annotation (Documentation(info="<html>

This function computes the Hessenberg matrix of matrix <b>A</b> by repetitive application of Householder similarity transformation
 <pre>
    <b>A</b>i+1 = (<b>I</b>-2*<b>u</b>_i*<b>u</b>_i')*<b>A</b>i*(<b>I</b>-2*<b>u</b>_i*<b>u</b>_i')
</pre>
with Householder vector <b>u</b>_i.
<p>
The elementary transformations can be subsumed under
 <pre> <b>A</b> -> <b>Q</b>*<b>A</b>*<b>Q</b>
</pre>
and <b>Q</b>*<b>A</b>*<b>Q</b> is Hessenberg matrix.
<p>
In contrast to function <b>Modelica_LinearSystems2.Math.Matrices.hess</b>, function <b>Modelica_LinearSystems2.Math.Matrices.hess3</b> does not use any LAPACK routine.

</html>"));
end hessenberg2;
