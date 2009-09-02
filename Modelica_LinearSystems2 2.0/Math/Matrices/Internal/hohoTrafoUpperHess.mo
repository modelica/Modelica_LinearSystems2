within Modelica_LinearSystems2.Math.Matrices.Internal;
function hohoTrafoUpperHess
  "Compute the similarity (Householder-) transformation S*A*S of matrix A with householder matrix S = I - 2u*u' to compute an upper Hessenberg form"

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

  Real P[na - r,na - r];
  Real A11[r,r]=A[1:r, 1:r];
  Real A12[r,na - r]=A[1:r, r + 1:na];
  Real A22[na - r,na - r]=A[r + 1:na, r + 1:na];

  Real alpha;

algorithm
  for i in 1:na loop
    S[i, i] := 1.0 + S[i, i];
                            //S=I-2u*u'
  end for;

  P := S[r + 1:na, r + 1:na];
  alpha := P[1, :]*A[r + 1:na, r];

  SAS := [A11,A12*P; [zeros(1, r - 1),matrix(alpha); zeros(na - r - 1, r)],P*
    A22*P];

  annotation (Documentation(info="<html>
<html>
This function calculates one step in the calculation of upper Hessenberg form. Therein it calculates the r'th column of the Hessenberg matrix which is of shape {x,x,fill(1,n-r-1)}.
</html>
</html>"));
end hohoTrafoUpperHess;
