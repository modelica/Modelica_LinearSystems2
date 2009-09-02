within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dgetrs
  "Solves a system of linear equations with the LU decomposition from dgetrf(..)"

  input Real LU[:,size(LU, 1)] "LU factorization of dgetrf of a square matrix";
  input Integer pivots[size(LU, 1)] "Pivot vector of dgetrf";
  input Real B[size(LU, 1),:] "Right hand side matrix B";
  output Real X[size(B, 1),size(B, 2)]=B "Solution matrix X";

  annotation (Documentation(info="
Purpose
=======
DGETRS solves a system of linear equations
   A * X = B  or  A' * X = B
with a general N-by-N matrix A using the LU factorization computed
by DGETRF.
Arguments
=========
TRANS   (input) CHARACTER*1
        Specifies the form of the system of equations:
        = 'N':  A * X = B  (No transpose)
        = 'T':  A'* X = B  (Transpose)
        = 'C':  A'* X = B  (Conjugate transpose = Transpose)
N       (input) INTEGER
        The order of the matrix A.  N >= 0.
NRHS    (input) INTEGER
        The number of right hand sides, i.e., the number of columns
        of the matrix B.  NRHS >= 0.
A       (input) DOUBLE PRECISION array, dimension (LDA,N)
        The factors L and U from the factorization A = P*L*U
        as computed by DGETRF.
LDA     (input) INTEGER
        The leading dimension of the array A.  LDA >= max(1,N).
IPIV    (input) INTEGER array, dimension (N)
        The pivot indices from DGETRF; for 1<=i<=N, row i of the
        matrix was interchanged with row IPIV(i).
B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
        On entry, the right hand side matrix B.
        On exit, the solution matrix X.
LDB     (input) INTEGER
        The leading dimension of the array B.  LDB >= max(1,N).
INFO    (output) INTEGER
        = 0:  successful exit
        < 0:  if INFO = -i, the i-th argument had an illegal value
"), Window(
      x=0.4,
      y=0.4,
      width=0.6,
      height=0.6));

protected
  Real work[size(LU, 1),size(LU, 1)]=LU;
  Integer info;
external "FORTRAN 77" dgetrs(
    "N",
    size(LU, 1),
    size(B, 2),
    work,
    size(LU, 1),
    pivots,
    X,
    size(B, 1),
    info)                 annotation (Library="lapack");
end dgetrs;
