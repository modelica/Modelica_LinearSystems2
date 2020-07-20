within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_singularValues
  "Compute singular values and left and right singular vectors of a complex matrix"
  import Complex;

  input Complex A[:,:] "Square or rectangular matrix";
  output Real sigma[min(size(A, 1), size(A, 2))] "singular values";
  output Complex U[size(A, 1),size(A, 1)] "Left orthogonal matrix";
  output Complex VT[size(A, 2),size(A, 2)] "Transposed right orthogonal matrix";

protected
  Integer info;
  Integer n=min(size(A, 1), size(A, 2)) "Number of singular values";

algorithm
  if n > 0 then
    (sigma,U,VT,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgesvd(
                                                                             A);
    assert(info == 0, "The numerical algorithm to compute the singular value decomposition did not converge");
  end if;
end C_singularValues;
