within Modelica_LinearSystems2.Math.Matrices;
function householderReflexion
  "reflect each of the vectors ai of matrix  A=[a1, a2, ..., an] on a plane with orthogonal vector u"
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[:,:];
  input Real u[size(A, 1)] "householder vector";

  output Real RA[size(A, 1),size(A, 2)] "reflexion of A";

protected
  Integer n=size(A, 2);
  Real h;
  Real lu=Vectors.length(u)*Vectors.length(u);

algorithm
  for i in 1:n loop
    h := scalar(2*transpose(matrix(u))*A[:, i]/lu);
    RA[:, i] := A[:, i] - h*u;
  end for;
end householderReflexion;
