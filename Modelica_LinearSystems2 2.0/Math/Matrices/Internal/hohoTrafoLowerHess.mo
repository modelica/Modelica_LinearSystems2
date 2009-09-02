within Modelica_LinearSystems2.Math.Matrices.Internal;
function hohoTrafoLowerHess
  "Compute the similarity transformation S*A*S of matrix A with householder matrix S = I - 2u*u' to compute a lower Hessenberg form"

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[:,size(A, 1)];
  input Real u[size(A, 1)] "householder vector";
  input Integer r;
  output Real SAS[size(A, 1),size(A, 1)];

protected
  Integer na=size(A, 1);
  Real S[:,:]=-2*matrix(u)*transpose(matrix(u))/(Vectors.length(u)*
      Vectors.length(u));   //S=u*u'/u'*u
  Integer i;

  Real P[na - r,na - r];
  Real A11[na - r,na - r]=A[1:na - r, 1:na - r];
  Real A21[r,na - r]=A[na - r + 1:na, 1:na - r];
  Real A22[r,r]=A[na - r + 1:na, na - r + 1:na];
  Real alpha;

algorithm
  for i in 1:na loop
    S[i, i] := 1.0 + S[i, i];   //S=I-2u*u'
  end for;

  P := S[1:na - r, 1:na - r];
  alpha := P[na - r, :]*A[1:na - r, na - r + 1];

  SAS := [P*A11*P,[zeros(na - r - 1, r); matrix(alpha),zeros(1, r - 1)]; A21*P,
    A22];

  annotation (Documentation(info="<html>
This function calculates one step in the calculation of lower Hessenberg form. Therein it calculates the (n-r+1)'th column of the Hessenberg matrix which is of shape {fill(1,n-r-1),x,x}.
</html>"));
end hohoTrafoLowerHess;
