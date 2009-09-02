within Modelica_LinearSystems2.Math.Matrices.Internal;
function householderSimilarityTransformation2
  "Calculate the similarity transformation SAS of matrix A with householder matrix S = I - 2u*u'/(u'*u) to compute a lower Hessenberg form"

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[:,size(A, 1)];
  input Real u[size(A, 1)] "householder vector";
  input Integer r;
  output Real SAS[size(A, 1),size(A, 1)];

protected
  Integer na=size(A, 1);
  Real S[:,:]=-2*matrix(u)*transpose(matrix(u))/(Vectors.length(u)*
      Vectors.length(u));                                                             //S=u*u'/u'*u
  Integer i;

algorithm
  assert(r <= na,
    "Input r in function \"Matrices.Internal.householderSimilarityTransformation2\" must fulfill r<=size(A,1)");
  for i in 1:na loop
    S[i, i] := 1.0 + S[i, i];   //S=I-2u*u'
  end for;

  SAS := if r < size(A, 1) then [S[1:r, 1:r]*A[1:r, 1:r]*S[1:r, 1:r],[
    zeros(r - 1, 1),A[1:r - 1, r + 2:na]; matrix(S[r, 1:r]*A[1:r, r + 1]),
    transpose(matrix(A[r, r + 2:na]))]; A[r + 1:na, 1:r]*S[1:r, 1:r],A[
    r + 1:na, r + 1:na]] else S*A*S;

  annotation (Documentation(info="<html>
It is assumed, that the input vector <b>u</b> is a Housholder vector of the shape
<blockquote><pre>
<b>u</b> = (u1, u2, ..., ur,0, ..., 0)
</pre></blockquote>
where r is an integer input. From
<blockquote><pre>
<b>S</b> = <b>I</b> - 2*<b>u</b>*<b>u</b>'/<b>u</b>'*<b>u</b> = [<b>P</b>, <b>0</b>; <b>0</b>, <b>I</b>]
</pre></blockquote>
with
<blockquote><pre>
dim(<b>P</b>) = r x r,   dim(<b>I</b>) = n-r x n-r
</pre></blockquote>
results
<blockquote><pre>
<b>S</b>*<b>A</b>*<b>S</b> = [<b>P</b>*<b>A</b>11*<b>P</b>, <b>P</b>*<b>A</b>12; <b>A</b>21*<b>P</b>, <b>A</b>22]
</pre></blockquote>
</html>"));
end householderSimilarityTransformation2;
