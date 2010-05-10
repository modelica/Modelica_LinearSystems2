within Modelica_LinearSystems2.Math.Matrices;
function householderSimilarityTransformation
  "Calculate the similarity transformation S*A*S of matrix A with symmetric householder matrix S = I - 2u*u'"

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[:,size(A, 1)] "Square matrix A";
  input Real u[size(A, 1)] "Householder vector";
  output Real SAS[size(A, 1),size(A, 1)];

protected
  Integer na=size(A, 1);
  Real S[:,:]=-2*matrix(u)*transpose(matrix(u))/(Vectors.length(u)*Vectors.length(u))
    "Symmetric matrix";
  Integer i;
algorithm
  for i in 1:na loop
    S[i, i] := 1.0 + S[i, i];
  end for;
  SAS := S*A*S;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<b>householderSimilarityTransformation</b>(A,u);
</pre></blockquote>
<h4>Description</h4>
<p>
This function computes the Housholder similarity transformation
</p>
<blockquote>
 <b>As</b> = <b>S</b>*<b>A</b>*<b>S</b>
</blockquote>
with
<blockquote>
 <b>S</b> = <b>I</b> -2*<b>u</b>*<b>u</b>'/(<b>u</b>'*<b>u</b>).
</blockquote>
<p>
This transformation is widely used for transforming non-symmetric matrices to a Hessenberg form.
</p>
<h4>Example</h4>
<blockquote><pre>
// First step of Hessenberg decomposition
  import   Modelica.Math.Vectors.Utilities;
  
  Real A[4,4] = [1,2,3,4;
                 3,4,5,6;
                 9,8,7,6;
                 1,2,0,0];
  Real Ar[4,4];
  Real u[4]={0,0,0,0};
  
  u[2:4]=Utilities.householderVector(A[2:4,1],{1,0,0});
  // u= = {0, 0.8107, 0.5819, 0.0647}
  Ar=householderSimilarityTransformation(A,u);
 //  Ar = [1.0,     -3.8787,    -1.2193,    3.531;
          -9.5394, 11.3407,      6.4336,   -5.9243;
           0.0,     3.1307,      0.7525,   -3.3670;
           0.0,     0.8021,     -1.1656,   -1.0932]
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.housholderReflection\">Matrices.housholderReflection</a>,
</p>
</html>"));
end householderSimilarityTransformation;
