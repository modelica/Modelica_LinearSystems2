within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dorglq
  "Generate a matrix Q with orthonormal rows which is defined as the product of elementary reflectors, as returned by DGELQF"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Orthogonal matrix of elementary reflectors";
  input Real tau[:] "Scalar factors of the elementary reflectors";

  output Real Q[size(A,1),size(A,2)]=A;
  output Integer info;
protected
  Integer m=size(A, 1);
  Integer n=size(A, 2);
  Integer k=size(tau, 1);
  Integer lda=max(1, m);
  Integer lwork=3*max(1, m);
  Real work[lwork];

external "Fortran 77" dorglq(
    m,
    n,
    k,
    Q,
    lda,
    tau,
    work,
    lwork,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
   which is defined as the first M rows of a product of K elementary
   reflectors of order N

         Q  =  H(k) . . . H(2) H(1)

   as returned by DGELQF.

   Arguments
   =========

   M       (input) INTEGER
           The number of rows of the matrix Q. M >= 0.

   N       (input) INTEGER
           The number of columns of the matrix Q. N >= M.

   K       (input) INTEGER
           The number of elementary reflectors whose product defines the
           matrix Q. M >= K >= 0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the i-th row must contain the vector which defines
           the elementary reflector H(i), for i = 1,2,...,k, as returned
           by DGELQF in the first k rows of its array argument A.
           On exit, the M-by-N matrix Q.

   LDA     (input) INTEGER
           The first dimension of the array A. LDA >= max(1,M).

   TAU     (input) DOUBLE PRECISION array, dimension (K)
           TAU(i) must contain the scalar factor of the elementary
           reflector H(i), as returned by DGELQF.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK. LWORK >= max(1,M).
           For optimum performance LWORK >= M*NB, where NB is
           the optimal blocksize.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument has an illegal value

   =====================================================================
"));
end dorglq;
