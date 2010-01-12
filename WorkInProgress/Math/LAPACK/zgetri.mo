within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zgetri
  import Modelica_LinearSystems2.Math.Complex;

  input Complex LU[:, size(LU, 1)]
    "LU factorization of zgetrf of a complex square matrix";
  input Integer pivots[size(LU, 1)] "Pivot vector of zgetrf";
  output Complex inv[size(LU, 1),size(LU, 2)];
  output Integer info;

protected
  Integer m=size(LU, 1);
  Integer n=size(LU, 2);
  Real LU_real[size(LU, 1),size(LU, 2)]=LU[:, :].re
    "LU factorization of zgetrf of a complex square matrix, real part";
  Real LU_imag[size(LU, 1),size(LU, 2)]=LU[:, :].im
    "LU factorization of zgetrf of a complex square matrix, imaginary part";
  Real inv_real[size(LU, 1),size(LU, 2)] "Inverse of matrix P*L*U, real part";
  Real inv_imag[size(LU, 1),size(LU, 2)]
    "Inverse of matrix P*L*U, imaginary part";

algorithm
  (inv_real,inv_imag,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zgetri(
                                                                                          LU_real, LU_imag, pivots);

  for l1 in 1:m loop
    for l2 in 1:n loop
      inv[l1, l2] := Complex(inv_real[l1, l2], inv_imag[l1, l2]);
    end for;
  end for;

end zgetri;
