within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zgesv "Solve complex system of linear equations A*X=B"
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;

  input Complex A[:,size(A, 1)];
  input Complex B[size(A, 1),:];
  output Complex X[size(A, 1),size(B, 2)]=B;
  output Integer info;

protected
  Integer n=size(A, 1) "Number of equations";

  Real A_real[size(A, 1),size(A, 2)]=A[:, :].re "Real part of matrix A";
  Real A_imag[size(A, 1),size(A, 2)]=A[:, :].im "Imaginary part of matrix A";
  Real B_real[size(B, 1),size(B, 2)]=B[:, :].re "Real part of matrix B";
  Real B_imag[size(B, 1),size(B, 2)]=B[:, :].im "Imaginary part of matrix B";
  Real X_real[size(X, 1),size(X, 2)] "Real part of matrix X";
  Real X_imag[size(X, 1),size(X, 2)] "Imaginary part of matrix X";

algorithm
  if n > 0 then
    (X_real,X_imag,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zgesv(
                                                          A_real, A_imag, B_real, B_imag);
    assert(info == 0, "Solving a linear system of equations with function
\"Matrices.C_solve2\" is not possible, because the system has either
no or infinitely many solutions (A is singular).");
    for l1 in 1:size(A, 1) loop
      for l2 in 1:size(B, 2) loop
        X[l1, l2] := Complex(X_real[l1, l2], X_imag[l1, l2]);
      end for;
    end for;
  end if;
end zgesv;
