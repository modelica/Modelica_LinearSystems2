within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_QR2
  "QR decomposition of a rectangular complex matrix without column pivoting (A = Q*R)"
  import Complex;

  input Complex A[:,:] "Rectangular matrix with size(A,1) >= size(A,2)";
  output Complex Q[size(A, 1),size(A, 1)]
    "Rectangular matrix with orthonormal columns such that Q*R=A[:,p]";
  output Complex R[size(A,1),size(A, 2)] "Square upper triangular matrix";

protected
  Complex A2[size(A, 1),size(A, 1)];
  Complex R2[size(A, 1),size(A, 1)];
  Integer nrow=size(A, 1);
  Integer ncol=size(A, 2);
  Integer minrowcol=min(nrow, ncol);

algorithm
  assert(nrow >= ncol, "\nInput matrix A[" + String(nrow) + "," + String(ncol) + "] has more columns as rows.
This is not allowed when calling Modelica.Matrices.C_QR(A).");

  if minrowcol > 0 then

    for l1 in 1:size(A, 1) loop
      for l2 in 1:size(A, 2) loop
        A2[l1, l2] := A[l1, l2];
      end for;
    end for;

    (Q,R2) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgeqrf(
                                                                  A2);
    for l1 in 1:size(R, 1) loop
      for l2 in l1:size(R, 2) loop
        R[l1, l2] := R2[l1, l2];
      end for;
    end for;

  else
    Q := fill(Complex(1), size(A, 1), 0);
    R := fill(Complex(0), 0, 0);
  end if;
end C_QR2;
