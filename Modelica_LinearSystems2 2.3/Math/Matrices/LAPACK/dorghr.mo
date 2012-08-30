within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dorghr
  "Generates a real orthogonal matrix Q which is defined as the product of IHI-ILO elementary reflectors of order N, as returned by DGEHRD"

  input Real A[:,size(A, 1)] "Square matrix A";
  input Integer ilo=1
    "lowest index where the original matrix had been Hessenbergform - ilo must have the same value as in the previous call of DGEHRD";
  input Integer ihi=size(A, 1)
    "highest index where the original matrix had been Hessenbergform  - ihi must have the same value as in the previous call of DGEHRD";
  input Real tau[max(0,size(A, 1) - 1)]
    "scalar factors of the elementary reflectors";
  output Real Aout[size(A, 1),size(A, 2)]=A
    "Orthogonal matrix as a result of elementary reflectors";
  output Integer info;
protected
  Integer n=size(A, 1);
  Integer lda=max(1, n);
  Integer lwork=max(1, 3*n);
  Real work[lwork];

external "Fortran 77" dorghr(
    n,
    ilo,
    ihi,
    Aout,
    lda,
    tau,
    work,
    lwork,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   DORGHR generates a real orthogonal matrix Q which is defined as the
   product of IHI-ILO elementary reflectors of order N, as returned by
   DGEHRD:

   Q = H(ilo) H(ilo+1) . . . H(ihi-1).

   Arguments
   =========

   N       (input) INTEGER
           The order of the matrix Q. N >= 0.

   ILO     (input) INTEGER
   IHI     (input) INTEGER
           ILO and IHI must have the same values as in the previous call
           of DGEHRD. Q is equal to the unit matrix except in the
           submatrix Q(ilo+1:ihi,ilo+1:ihi).
           1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
           On entry, the vectors which define the elementary reflectors,
           as returned by DGEHRD.
           On exit, the N-by-N orthogonal matrix Q.

   LDA     (input) INTEGER
           The leading dimension of the array A. LDA >= max(1,N).

   TAU     (input) DOUBLE PRECISION array, dimension (N-1)
           TAU(i) must contain the scalar factor of the elementary
           reflector H(i), as returned by DGEHRD.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK. LWORK >= IHI-ILO.
           For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
           the optimal blocksize.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   =====================================================================
"));
end dorghr;
