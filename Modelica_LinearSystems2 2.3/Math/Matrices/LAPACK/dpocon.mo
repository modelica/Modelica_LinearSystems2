within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dpocon
  "Estimates the reciprocal of the condition number (1-norm) of a real symmetric matrix A using the Cholesky factor"

  input Real cholA[:,size(cholA,1)] "Cholesky factor of matrix A";
  input Real anorm "Norm of A";
  input Boolean upper "Specifies whether A is upper or lower Cholesky factor";
  output Real rcond "Reciprocal condition number of A";
  output Integer info;
protected
  String uplo=if upper then "U" else "L";
  Integer n=size(cholA, 1);
  Integer lda=max(1,n);
  Real work[3*n];
  Integer iwork[n];

external "Fortran 77" dpocon(
    uplo,
    n,
    cholA,
    lda,
    anorm,
    rcond,
    work,
    iwork,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   DPOCON estimates the reciprocal of the condition number (in the
   1-norm) of a real symmetric positive definite matrix using the
   Cholesky factorization A = U**T*U or A = L*L**T computed by DPOTRF.

   An estimate is obtained for norm(inv(A)), and the reciprocal of the
   condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

   Arguments
   =========

   UPLO    (input) CHARACTER*1
           = 'U':  Upper triangle of A is stored;
           = 'L':  Lower triangle of A is stored.

   N       (input) INTEGER
           The order of the matrix A.  N >= 0.

   A       (input) DOUBLE PRECISION array, dimension (LDA,N)
           The triangular factor U or L from the Cholesky factorization
           A = U**T*U or A = L*L**T, as computed by DPOTRF.

   LDA     (input) INTEGER
           The leading dimension of the array A.  LDA >= max(1,N).

   ANORM   (input) DOUBLE PRECISION
           The 1-norm (or infinity-norm) of the symmetric matrix A.

   RCOND   (output) DOUBLE PRECISION
           The reciprocal of the condition number of the matrix A,
           computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
           estimate of the 1-norm of inv(A) computed in this routine.

   WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)

   IWORK   (workspace) INTEGER array, dimension (N)

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   =====================================================================


"));
end dpocon;
