within Modelica_LinearSystems2.Math.Matrices.Internal;
function householderSimilarityTransformation2
  "Calculate the similarity transformation SAS of matrix A with householder matrix S = I - 2u*u'/(u'*u) to compute a lower Hessenberg form"

  import Modelica.Math.Vectors.length;

  input Real A[:,size(A, 1)];
  input Real u[size(A, 1)] "householder vector";
  input Integer r;
  output Real SAS[size(A, 1),size(A, 1)];

protected
  Integer na=size(A, 1);
  Real S[:,:]=-2*matrix(u)*transpose(matrix(u))/(length(u)*length(u)); //S=u*u'/u'*u
  Integer i;

algorithm
  assert(r <= na,
    "Input r in function \"Modelica_LinearSystems2.Math.Matrices.Internal.householderSimilarityTransformation2\" must fulfill r<=size(A,1)");
  for i in 1:na loop
    S[i, i] := 1.0 + S[i, i];   //S=I-2u*u'
  end for;

  SAS := if r < size(A, 1) then [S[1:r, 1:r]*A[1:r, 1:r]*S[1:r, 1:r],[
    zeros(r - 1, 1),A[1:r - 1, r + 2:na]; matrix(S[r, 1:r]*A[1:r, r + 1]),
    transpose(matrix(A[r, r + 2:na]))]; A[r + 1:na, 1:r]*S[1:r, 1:r],A[
    r + 1:na, r + 1:na]] else S*A*S;

  annotation (Documentation(info="<html>
<p>
It is assumed that the input vector <strong>u</strong> is a Housholder vector of the shape
</p>
<blockquote><pre>
<strong>u</strong> = (u1, u2, ..., ur,0, ..., 0)
</pre></blockquote>
<p>
where r is an integer input. From
</p>
<blockquote><pre>
<strong>S</strong> = <strong>I</strong> - 2*<strong>u</strong>*<strong>u</strong>'/<strong>u</strong>'*<strong>u</strong> = [<strong>P</strong>, <strong>0</strong>; <strong>0</strong>, <strong>I</strong>]
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
dim(<strong>P</strong>) = r x r,   dim(<strong>I</strong>) = n-r x n-r
</pre></blockquote>
<p>
results
</p>
<blockquote><pre>
<strong>S</strong>*<strong>A</strong>*<strong>S</strong> = [<strong>P</strong>*<strong>A</strong>11*<strong>P</strong>, <strong>P</strong>*<strong>A</strong>12; <strong>A</strong>21*<strong>P</strong>, <strong>A</strong>22]
</pre></blockquote>
</html>"));
end householderSimilarityTransformation2;
