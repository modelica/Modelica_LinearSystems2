within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zgesvd "Determine singular values of a complex matrix"
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;

  input Complex A[:,:] "Square or rectangular matrix";
  output Real sigma[min(size(A, 1), size(A, 2))] "singular values";
  output Complex U[size(A, 1),size(A, 1)] "Left orthogonal matrix";
  output Complex VT[size(A, 2),size(A, 2)] "Transposed right orthogonal matrix";
  output Integer info;

protected
  Integer n=min(size(A, 1), size(A, 2)) "Number of singular values";

  Real U_real[size(A, 1),size(A, 1)] "Left orthogonal matrix, real part";
  Real U_imag[size(A, 1),size(A, 1)] "Left orthogonal matrix, imaginary part";
  Real VT_real[size(A, 2),size(A, 2)] "Left orthogonal matrix, real part";
  Real VT_imag[size(A, 2),size(A, 2)] "Left orthogonal matrix, imaginary part";
  Real A_real[size(A, 1),size(A, 2)]=A[:, :].re "Real part of matrix A";
  Real A_imag[size(A, 1),size(A, 2)]=A[:, :].im "Imaginary part of matrix A";

algorithm
  if n > 0 then

    (sigma,U_real,U_imag,VT_real,VT_imag,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zgesvd(
                                                                                 A_real, A_imag);
    assert(info == 0, "The numerical algorithm to compute the singular value decomposition did not converge");

    for l1 in 1:size(A, 1) loop
      for l2 in 1:size(A, 1) loop
        U[l1, l2] := Complex(U_real[l1, l2], U_imag[l1, l2]);
      end for;
    end for;
    for l1 in 1:size(A, 2) loop
      for l2 in 1:size(A, 2) loop
        VT[l1, l2] := Complex(VT_real[l1, l2], VT_imag[l1, l2]);
      end for;
    end for;

  end if;
end zgesvd;
