within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dggev_eigenValues
  "Compute only generalized eigenvalues for a (A,B) system"
  extends Modelica.Icons.Function;

  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),size(A, 1)];
  output Real alphaReal[size(A, 1)]
    "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag)/beta)";
  output Real alphaImag[size(A, 1)] "Imaginary part of alpha";
  output Real beta[size(A, 1)] "Denominator of eigenvalue";
  output Integer info;
protected
  Integer n=size(A, 1);
  Integer lda=max(1,n);
  Integer lwork=max(1,2*8*n);
  Real Awork[n,n]=A;
  Real Bwork[n,n]=B;
  Real work[lwork];
  Real lEigenVectors[1];
  Real rEigenVectors[1];

external "Fortran 77" dggev(
    "N",
    "N",
    size(A,1),
    Awork,
    lda,
    Bwork,
    lda,
    alphaReal,
    alphaImag,
    beta,
    lEigenVectors,
    1,
    rEigenVectors,
    1,
    work,
    lwork,
    info) annotation(Library = {"lapack"});
  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   the generalized eigenvalues, and optionally, the left and/or right
   generalized eigenvectors.

   A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   singular. It is usually represented as the pair (alpha,beta), as
   there is a reasonable interpretation for beta=0, and even for both
   being zero.

   The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   of (A,B) satisfies

                    A * v(j) = lambda(j) * B * v(j).

   The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   of (A,B) satisfies

                    u(j)**H * A  = lambda(j) * u(j)**H * B .

   where u(j)**H is the conjugate-transpose of u(j).

   Arguments
   =========

   JOBVL   (input) CHARACTER*1
           = 'N':  do not compute the left generalized eigenvectors;
           = 'V':  compute the left generalized eigenvectors.

   JOBVR   (input) CHARACTER*1
           = 'N':  do not compute the right generalized eigenvectors;
           = 'V':  compute the right generalized eigenvectors.

   N       (input) INTEGER
           The order of the matrices A, B, VL, and VR.  N >= 0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
           On entry, the matrix A in the pair (A,B).
           On exit, A has been overwritten.

   LDA     (input) INTEGER
           The leading dimension of A.  LDA >= max(1,N).

   B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
           On entry, the matrix B in the pair (A,B).
           On exit, B has been overwritten.

   LDB     (input) INTEGER
           The leading dimension of B.  LDB >= max(1,N).

   ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   BETA    (output) DOUBLE PRECISION array, dimension (N)
           On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
           be the generalized eigenvalues.  If ALPHAI(j) is zero, then
           the j-th eigenvalue is real; if positive, then the j-th and
           (j+1)-st eigenvalues are a complex conjugate pair, with
           ALPHAI(j+1) negative.

           Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
           may easily over- or underflow, and BETA(j) may even be zero.
           Thus, the user should avoid naively computing the ratio
           alpha/beta.  However, ALPHAR and ALPHAI will be always less
           than and usually comparable with norm(A) in magnitude, and
           BETA always less than and usually comparable with norm(B).

   VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
           If JOBVL = 'V', the left eigenvectors u(j) are stored one
           after another in the columns of VL, in the same order as
           their eigenvalues. If the j-th eigenvalue is real, then
           u(j) = VL(:,j), the j-th column of VL. If the j-th and
           (j+1)-th eigenvalues form a complex conjugate pair, then
           u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
           Each eigenvector will be scaled so the largest component have
           abs(real part)+abs(imag. part)=1.
           Not referenced if JOBVL = 'N'.

   LDVL    (input) INTEGER
           The leading dimension of the matrix VL. LDVL >= 1, and
           if JOBVL = 'V', LDVL >= N.

   VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
           If JOBVR = 'V', the right eigenvectors v(j) are stored one
           after another in the columns of VR, in the same order as
           their eigenvalues. If the j-th eigenvalue is real, then
           v(j) = VR(:,j), the j-th column of VR. If the j-th and
           (j+1)-th eigenvalues form a complex conjugate pair, then
           v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
           Each eigenvector will be scaled so the largest component have
           abs(real part)+abs(imag. part)=1.
           Not referenced if JOBVR = 'N'.

   LDVR    (input) INTEGER
           The leading dimension of the matrix VR. LDVR >= 1, and
           if JOBVR = 'V', LDVR >= N.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK.  LWORK >= max(1,8*N).
           For good performance, LWORK must generally be larger.

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value.
           = 1,...,N:
                 The QZ iteration failed.  No eigenvectors have been
                 calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
                 should be correct for j=INFO+1,...,N.
           > N:  =N+1: other than QZ iteration failed in DHGEQZ.
                 =N+2: error return from DTGEVC.

   =====================================================================
"));
end dggev_eigenValues;
