within Modelica_LinearSystems2.Math.Matrices.Internal;
function reorderRSF2
  "Reorders a real Schur factorization for poleAssignmentMI design"

  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real T[:,:] "upper quasi-triangular matrix in Schur canonical form";
  input Real Q[:,size(T, 2)] "matrix of Schur vectors";
  input Real alphaReal[size(T, 1)]
    "Real part of eigenvalue=alphaReal+i*alphaImag";
  input Real alphaImag[size(T, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";
  input Real alpha
    "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";

  output Real To[size(T, 1),size(T, 2)];
  output Real Qo[size(T, 1),size(T, 2)];
  output Real wr[size(T, 2)];
  output Real wi[size(T, 2)];

protected
  Integer n=size(T, 2);
  Boolean select[n]=fill(false, n);
  Integer i;
algorithm
  for i in 1:n loop
    if alphaReal[i] < alpha then
      select[i] := true;
    end if;
  end for;

  (To,Qo,wr,wi) := LAPACK.dtrsen(
      "E",
      "V",
      select,
      T,
      Q);

end reorderRSF2;
