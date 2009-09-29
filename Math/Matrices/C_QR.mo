within Modelica_LinearSystems2.Math.Matrices;
function C_QR
  "QR decomposition of a rectangular complex matrix without column pivoting (A = Q*R)"
  import Modelica_LinearSystems2.Math.Complex;

  input Complex A[:,:] "Rectangular matrix with size(A,1) >= size(A,2)";
  output Complex Q[size(A, 1),size(A, 2)]
    "Rectangular matrix with orthonormal columns such that Q*R=A[:,p]";
  output Complex R[min(size(A, 1), size(A, 2)),size(A, 2)]
    "Square upper triangular matrix";

protected
  Integer nrow=size(A, 1);
  Integer ncol=size(A, 2);
  Integer minrowcol=min(nrow, ncol);

algorithm
  assert(nrow >= ncol, "\nInput matrix A[" + String(nrow) + "," + String(ncol)
     +                                                                     "] has more columns as rows.
This is not allowed when calling Modelica.Matrices.C_QR(A).");

  if minrowcol > 0 then

    (Q,R) := Modelica_LinearSystems2.Math.Matrices.LAPACK.zgeqrf(A);

  else
    Q := fill(
      Complex(1),
      size(A, 1),
      0);
    R := fill(
      Complex(0),
      0,
      0);
  end if;
end C_QR;
