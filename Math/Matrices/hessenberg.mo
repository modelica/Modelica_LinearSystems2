within Modelica_LinearSystems2.Math.Matrices;
function hessenberg
  "Compute an upper Hessenberg matrix by repeatedly applicated householder similarity transformation"
  import Modelica;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real H[size(H, 1),size(H, 2)];

  output Real Ht[size(H, 1),size(H, 2)];

protected
  Integer q=size(H, 1);
  Real u[q] "householder vector";
  Integer ll;

algorithm
  Ht := Modelica_LinearSystems2.Math.Matrices.toUpperHessenberg(
    H,
    1,
    size(H, 1));

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
<b>Modelica_LinearSystems2.Math.Matrices.hess</b> uses LAPACK routine dgehrd. In contrast to this function <b>Modelica_LinearSystems2.Math.Matrices.Internal.hessenberg</b> does not use any LAPACK routine.
 
</html>"));
end hessenberg;
