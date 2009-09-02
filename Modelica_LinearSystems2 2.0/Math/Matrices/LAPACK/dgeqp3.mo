within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dgeqp3 "computes a QR factorization with column pivoting"

  input Real A[:,:];
  input Integer lwork1=3*size(A, 2) + 1
    "size of work array; should be optimized with Modelica_LinearSystems2.Math.Matrices.Internal.dgeqp3_workdim";
  output Real Aout[size(A, 1),size(A, 2)]=A
    "the upper triangle of the array contains the upper trapezoidal matrix R; the elements below the diagonal, together with the array TAU, represent the orthogonal matrix Q as a product of elementary reflectors";
  output Integer jpvt[size(A, 2)] "pivoting indices";
  output Real tau[min(size(A, 1), size(A, 2))]
    "scalar factors of the elementary reflectors";
  output Integer info;
  output Real work[max(lwork1, 3*size(A, 2) + 1)];
protected
  Integer m=size(A, 1);
  Integer n=size(A, 2);
  Integer lda=max(1, m);
  Integer lwork2=if lwork1 == -1 then -1 else max(1, lwork1);

  annotation (Documentation(info="   Purpose
   =======

   DGEQP3 computes a QR factorization with column pivoting of a
   matrix A:  A*P = Q*R  using Level 3 BLAS.

   Arguments
   =========

   M       (input) INTEGER
           The number of rows of the matrix A. M >= 0.

   N       (input) INTEGER
           The number of columns of the matrix A.  N >= 0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the M-by-N matrix A.
           On exit, the upper triangle of the array contains the
           min(M,N)-by-N upper trapezoidal matrix R; the elements below
           the diagonal, together with the array TAU, represent the
           orthogonal matrix Q as a product of min(M,N) elementary
           reflectors.

   LDA     (input) INTEGER
           The leading dimension of the array A. LDA >= max(1,M).

   JPVT    (input/output) INTEGER array, dimension (N)
           On entry, if JPVT(J).ne.0, the J-th column of A is permuted
           to the front of A*P (a leading column); if JPVT(J)=0,
           the J-th column of A is a free column.
           On exit, if JPVT(J)=K, then the J-th column of A*P was the
           the K-th column of A.

   TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
           The scalar factors of the elementary reflectors.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
           On exit, if INFO=0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK. LWORK >= 3*N+1.
           For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
           is the optimal blocksize.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0: successful exit.
           < 0: if INFO = -i, the i-th argument had an illegal value.

   Further Details
   ===============

   The matrix Q is represented as a product of elementary reflectors

      Q = H(1) H(2) . . . H(k), where k = min(m,n).

   Each H(i) has the form

      H(i) = I - tau * v * v'

   where tau is a real/complex scalar, and v is a real/complex vector
   with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
   A(i+1:m,i), and tau in TAU(i).

   Based on contributions by
     G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
     X. Sun, Computer Science Dept., Duke University, USA

   =====================================================================
"));
external "Fortran 77" dgeqp3(
    m,
    n,
    Aout,
    lda,
    jpvt,
    tau,
    work,
    lwork2,
    info) annotation(Library = {"lapack"});

end dgeqp3;
