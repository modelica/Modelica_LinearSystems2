within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zgetrf
  import Complex;

  input Complex A[:,:] "Square or rectangular matrix";
  output Complex LU[size(A, 1),size(A, 2)];
  output Integer pivots[min(size(A, 1), size(A, 2))] "Pivot vector";
  output Integer info "Information";

protected
  Integer m=size(A, 1);
  Integer n=size(A, 2);
  Real A_real[size(A, 1),size(A, 2)]=A[:, :].re "Square or rectangular matrix";
  Real A_imag[size(A, 1),size(A, 2)]=A[:, :].im "Square or rectangular matrix";
  Real LU_real[size(A, 1),size(A, 2)]
    "LU factorization in packed format, real part";
  Real LU_imag[size(A, 1),size(A, 2)]
    "LU factorization in packed format, imaginary part";

algorithm
  (LU_real,LU_imag,pivots,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zgetrf(
                                                                                            A_real, A_imag);
  for l1 in 1:m loop
    for l2 in 1:n loop
      LU[l1, l2] := Complex(LU_real[l1, l2], LU_imag[l1, l2]);
    end for;
  end for;

end zgetrf;
