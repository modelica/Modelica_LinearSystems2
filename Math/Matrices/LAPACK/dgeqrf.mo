within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dgeqrf "Computes a QR factorization without pivoting"

  input Real A[:,:];
  input Integer lwork1=size(A, 2)
    "size of work array; should be optimized with Modelica_LinearSystems2.Math.Matrices.Internal.dgeqp3_workdim";
  output Real Aout[size(A, 1),size(A, 2)]=A
    "the upper triangle of the array contains the upper trapezoidal matrix R; the elements below the diagonal, together with the array TAU, represent the orthogonal matrix Q as a product of elementary reflectors";
  output Real tau[min(size(A, 1), size(A, 2))]
    "scalar factors of the elementary reflectors";
  output Integer info;
  output Real work[max(lwork1, 3*size(A, 2) + 1)];
protected
  Integer m=size(A, 1);
  Integer n=size(A, 2);
  Integer lda=max(1, m);
  Integer lwork2=if lwork1 == -1 then -1 else max(1, lwork1);

external "Fortran 77" dgeqrf(
    m,
    n,
    Aout,
    lda,
    tau,
    work,
    lwork2,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="
   Purpose
   =======

   DGEQRF computes a QR factorization of a real M-by-N matrix A:
   A = Q * R.

   Arguments
   =========

   M       (input) INTEGER
           The number of rows of the matrix A.  M >= 0.

   N       (input) INTEGER
           The number of columns of the matrix A.  N >= 0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the M-by-N matrix A.
           On exit, the elements on and above the diagonal of the array
           contain the min(M,N)-by-N upper trapezoidal matrix R (R is
           upper triangular if m >= n); the elements below the diagonal,
           with the array TAU, represent the orthogonal matrix Q as a
           product of min(m,n) elementary reflectors (see Further
           Details).

   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,M).

   TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
           The scalar factors of the elementary reflectors (see Further
           Details).

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK.  LWORK >= max(1,N).
           For optimum performance LWORK >= N*NB, where NB is
           the optimal blocksize.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   Further Details
   ===============

   The matrix Q is represented as a product of elementary reflectors

      Q = H(1) H(2) . . . H(k), where k = min(m,n).

   Each H(i) has the form

      H(i) = I - tau * v * v'

   where tau is a real scalar, and v is a real vector with
   v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
   and tau in TAU(i).

   =====================================================================  "));
end dgeqrf;
