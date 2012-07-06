within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function dtrevc
  "Compute the right and/or left eigenvectors of a real upper quasi-triangular matrix T"
  input Real T[:,size(T, 1)];
  input String side="R";
  input String howmny="B";
  input Real Q[size(T, 1),size(T, 1)];

  output Real lEigenVectors[size(T, 1),size(T, 1)]=Q
    "left eigenvectors of matrix T";
  output Real rEigenVectors[size(T, 1),size(T, 1)]=Q
    "right eigenvectors of matrix T";
  output Integer info;

protected
  Integer n=size(T, 1);
  Boolean select[n];
  Integer ldt=max(1, n);
  Integer ldvl=max(1, n);
  Integer ldvr=max(1, n);
  Real work[3*n];

external "Fortran 77" dtrevc(
    side,
    howmny,
    select,
    n,
    T,
    ldt,
    lEigenVectors,
    ldvl,
    rEigenVectors,
    ldvr,
    n,
    n,
    work,
    info) annotation(Library = {"lapack"});

  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   DTREVC computes some or all of the right and/or left eigenvectors of
   a real upper quasi-triangular matrix T.
   Matrices of this type are produced by the Schur factorization of
   a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.

   The right eigenvector x and the left eigenvector y of T corresponding
   to an eigenvalue w are defined by:

      T*x = w*x,     (y**H)*T = w*(y**H)

   where y**H denotes the conjugate transpose of y.
   The eigenvalues are not input to this routine, but are read directly
   from the diagonal blocks of T.

   This routine returns the matrices X and/or Y of right and left
   eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
   input matrix.  If Q is the orthogonal factor that reduces a matrix
   A to Schur form T, then Q*X and Q*Y are the matrices of right and
   left eigenvectors of A.

   Arguments
   =========

   SIDE    (input) CHARACTER*1
           = 'R':  compute right eigenvectors only;
           = 'L':  compute left eigenvectors only;
           = 'B':  compute both right and left eigenvectors.

   HOWMNY  (input) CHARACTER*1
           = 'A':  compute all right and/or left eigenvectors;
           = 'B':  compute all right and/or left eigenvectors,
                   backtransformed by the matrices in VR and/or VL;
           = 'S':  compute selected right and/or left eigenvectors,
                   as indicated by the logical array SELECT.

   SELECT  (input/output) LOGICAL array, dimension (N)
           If HOWMNY = 'S', SELECT specifies the eigenvectors to be
           computed.
           If w(j) is a real eigenvalue, the corresponding real
           eigenvector is computed if SELECT(j) is .TRUE..
           If w(j) and w(j+1) are the real and imaginary parts of a
           complex eigenvalue, the corresponding complex eigenvector is
           computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
           on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
           .FALSE..
           Not referenced if HOWMNY = 'A' or 'B'.

   N       (input) INTEGER
           The order of the matrix T. N >= 0.

   T       (input) DOUBLE PRECISION array, dimension (LDT,N)
           The upper quasi-triangular matrix T in Schur canonical form.

   LDT     (input) INTEGER
           The leading dimension of the array T. LDT >= max(1,N).

   VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
           On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
           contain an N-by-N matrix Q (usually the orthogonal matrix Q
           of Schur vectors returned by DHSEQR).
           On exit, if SIDE = 'L' or 'B', VL contains:
           if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
           if HOWMNY = 'B', the matrix Q*Y;
           if HOWMNY = 'S', the left eigenvectors of T specified by
                            SELECT, stored consecutively in the columns
                            of VL, in the same order as their
                            eigenvalues.
           A complex eigenvector corresponding to a complex eigenvalue
           is stored in two consecutive columns, the first holding the
           real part, and the second the imaginary part.
           Not referenced if SIDE = 'R'.

   LDVL    (input) INTEGER
           The leading dimension of the array VL.  LDVL >= 1, and if
           SIDE = 'L' or 'B', LDVL >= N.

   VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
           On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
           contain an N-by-N matrix Q (usually the orthogonal matrix Q
           of Schur vectors returned by DHSEQR).
           On exit, if SIDE = 'R' or 'B', VR contains:
           if HOWMNY = 'A', the matrix X of right eigenvectors of T;
           if HOWMNY = 'B', the matrix Q*X;
           if HOWMNY = 'S', the right eigenvectors of T specified by
                            SELECT, stored consecutively in the columns
                            of VR, in the same order as their
                            eigenvalues.
           A complex eigenvector corresponding to a complex eigenvalue
           is stored in two consecutive columns, the first holding the
           real part and the second the imaginary part.
           Not referenced if SIDE = 'L'.

   LDVR    (input) INTEGER
           The leading dimension of the array VR.  LDVR >= 1, and if
           SIDE = 'R' or 'B', LDVR >= N.

   MM      (input) INTEGER
           The number of columns in the arrays VL and/or VR. MM >= M.

   M       (output) INTEGER
           The number of columns in the arrays VL and/or VR actually
           used to store the eigenvectors.
           If HOWMNY = 'A' or 'B', M is set to N.
           Each selected real eigenvector occupies one column and each
           selected complex eigenvector occupies two columns.

   WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)

   INFO    (output) INTEGER
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value

   Further Details
   ===============

   The algorithm used in this program is basically backward (forward)
   substitution, with scaling to make the the code robust against
   possible overflow.

   Each eigenvector is normalized so that the element of largest
   magnitude has magnitude 1; here the magnitude of a complex number
   (x,y) is taken to be |x| + |y|.

   =====================================================================  "));
end dtrevc;
