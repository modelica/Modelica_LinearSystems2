within Modelica_LinearSystems2.Math.Matrices;
function householderReflexion
  "Reflect each of the vectors ai of matrix  A=[a1, a2, ..., an] on a plane with orthogonal vector u"
  import Modelica.Math.Vectors.length;

  input Real A[:,:] "Rectangular matrix";
  input Real u[size(A, 1)] "Householder vector";

  output Real RA[size(A, 1),size(A, 2)] "Reflexion of A";

protected
  Integer n=size(A, 2);
  Real h;
  Real lu=length(u)*length(u);

algorithm
  for i in 1:n loop
    h := scalar(2*transpose(matrix(u))*A[:, i]/lu);
    RA[:, i] := A[:, i] - h*u;
  end for;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<b>householderReflection</b>(A,u);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the Housholder reflection (transformation)
</p>
<blockquote>
<b>Ar</b> = <b>Q</b>*<b>A</b>
</blockquote>
<p>
with
</p>
<blockquote>
<b>Q</b> = <b>I</b> -2*<b>u</b>*<b>u</b>'/(<b>u</b>'*<b>u</b>)
</blockquote>
<p>
where <b>u</b>*<b>u</b> is housholder vector, i.e. the normal vector of the reflection plane.
</p>
<p>
Householder reflection is widely used in numerical linear algebra, e.g. to perform QR decompositions.
</p>

<h4>Example</h4>
<blockquote><pre>
// First step of QR decomposition
  import   Modelica.Math.Vectors.Utilities;

  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real Ar[3,3];
  Real u[:];

  u=Utilities.householderVector(A[:,1],{1,0,0});
  // u={ = {0.763, 0.646, 0}
  Ar=householderReflexion(A,u);
 // Ar = [-6.0828,   -5.2608,   -4.4388;
 //        0.0,      -1.1508,   -2.3016;
 //        0.0,       2.0,       0.0]

</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.Utilities.householderSimilarityTransformation\">Matrices.housholderSimilarityTransformation</a>
</p>
</html>"));
end householderReflexion;
