within ;
package ObsoleteLinearSystems2
  "Library that contains components from Modelica_LinearSystems2 Library 2.4.X that have been removed from version 3.0.0"

  package Math "Package of additional functions for Modelica.Math"
    extends Modelica.Icons.Package;

    operator record Complex "Record defining a Complex number"

      encapsulated function j "Obsolete imaginary unit function"
        import Modelica;
        import Complex;

        output Complex c = Modelica.ComplexMath.j "= sqrt(-1)";

        annotation(
          obsolete = "Obsolete function. See function documentation for how to migrate imaginary unit j.",
          Inline=true,
          Documentation(info="<html>
<p>
A&nbsp;definition of complex number <code>c1</code> such as
</p>
<blockquote><pre>
Complex j = Modelica_LinearSystems2.Math.Complex.j();
Complex c1=2+3*j;
</pre></blockquote>
<p>
using Modelica_LinearSystems2 <strong>2.4.X</strong> shall be
automatically migrated now to (just the first line is affected)
</p>
<blockquote><pre>
Complex j = ObsoleteLinearSystems2.Math.Complex.j();
Complex c1=2+3*j;
</pre></blockquote>

<p>
The automatic conversion to proper <code>Modelica.ComplexMath.j</code> is
not possible since this is a&nbsp;complex constant. Therefore, you shall replace
the first line manually to
</p>
<blockquote><pre>
Complex j = Modelica.ComplexMath.j;
Complex c1=2+3*j;  // this stays unchanged
</pre></blockquote>
<p>
if <code>j</code> is further needed in the class.
In many classes, <code>j</code> is used only for the abovementioned complex
number definion. Then, the import of <code>j</code> can be completely omitted and
only the second line can be modified as follows
</p>
<blockquote><pre>
Complex c1 = Complex(2, 3);
</pre></blockquote>
</html>"));
      end j;

    end Complex;

    package LAPACK "Package of LAPACK functions"
      extends Modelica.Icons.Package;

      function dgeev
        "Obsolete: Compute the eigenvalues and the (real) left and right eigenvectors of matrix A"

        input Real A[:,size(A, 1)];

        output Real alphaReal[size(A, 1)]
          "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
        output Real alphaImag[size(A, 1)]
          "Imaginary part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
        output Real lEigenVectors[size(A, 1),size(A, 1)]
          "left eigenvectors of matrix A";
        output Real rEigenVectors[size(A, 1),size(A, 1)]
          "right eigenvectors of matrix A";
        output Integer info;
      protected
        Integer n=size(A, 1);
        Integer lwork=4*n;
        Real work[lwork];

      external "Fortran 77" dgeev(
          "V",
          "V",
          n,
          A,
          n,
          alphaReal,
          alphaImag,
          lEigenVectors,
          n,
          rEigenVectors,
          n,
          work,
          size(work, 1),
          info) annotation(Library = {"lapack"});

        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dgeev instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dgeev\">Modelica.Math.Matrices.LAPACK.dgeev</a>
instead.
</p>
<p>
Note: output <code>lEigenVector</code> is missing in the function from the Modelica
Standard Library and must be removed or resolved in another way.
</p>
</html>"));
      end dgeev;

      function dgegv "Compute generalized eigenvalues for a (A,B) system"

        input Real A[:,size(A, 1)];
        input Real B[size(A, 1),size(A, 1)];
        output Real alphaReal[size(A, 1)]
          "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag)/beta)";
        output Real alphaImag[size(A, 1)] "Imaginary part of alpha";
        output Real beta[size(A, 1)] "Denominator of eigenvalue";
        output Integer info;
      protected
        Integer n=size(A, 1);
        Integer lwork=2*n + max(6*n, n*n + n);
        Real Awork[n,n]=A;
        Real Bwork[n,n]=B;
        Real work[lwork];
        Real dummy1[1,1];
        Real dummy2[1,1];

        external "Fortran 77" dgegv(
          "N",
          "N",
          n,
          Awork,
          n,
          Bwork,
          n,
          alphaReal,
          alphaImag,
          beta,
          dummy1,
          1,
          dummy2,
          1,
          work,
          size(work, 1),
          info) annotation (Library="lapack");
        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dggev instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dggev\">Modelica.Math.Matrices.LAPACK.dggev</a>
instead, see also below the documentation of the original routine DGEGV.
</p>
<p>
Note: there are two additional outputs in the abovementioned function
from the Modelica Standard Library: <code>lEigenVectors</code> and
<code>rEigenVectors</code>.
They must be added in a&nbsp;function call or resolved in another way.
</p>

<pre>
Lapack documentation:

   Purpose
   =======

   This routine is deprecated and has been replaced by routine DGGEV.

   DGEGV computes for a pair of n-by-n real nonsymmetric matrices A and
   B, the generalized eigenvalues (alphar +/- alphai*i, beta), and
   optionally, the left and/or right generalized eigenvectors (VL and
   VR).

   A generalized eigenvalue for a pair of matrices (A,B) is, roughly
   speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
   is singular.  It is usually represented as the pair (alpha,beta),
   as there is a reasonable interpretation for beta=0, and even for
   both being zero.  A good beginning reference is the book, \"Matrix
   Computations\", by G. Golub &amp; C. van Loan (Johns Hopkins U. Press)

   A right generalized eigenvector corresponding to a generalized
   eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such
   that  (A - w B) r = 0 .  A left generalized eigenvector is a vector
   l such that l**H * (A - w B) = 0, where l**H is the
   conjugate-transpose of l.

   Note: this routine performs \"full balancing\" on A and B -- see
   \"Further Details\", below.

   Arguments
   =========

   JOBVL   (input) CHARACTER*1
           = 'N':  do not compute the left generalized eigenvectors;
           = 'V':  compute the left generalized eigenvectors.

   JOBVR   (input) CHARACTER*1
           = 'N':  do not compute the right generalized eigenvectors;
           = 'V':  compute the right generalized eigenvectors.

   N       (input) INTEGER
           The order of the matrices A, B, VL, and VR.  N &gt;= 0.

   A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
           On entry, the first of the pair of matrices whose
           generalized eigenvalues and (optionally) generalized
           eigenvectors are to be computed.
           On exit, the contents will have been destroyed.  (For a
           description of the contents of A on exit, see \"Further
           Details\", below.)

   LDA     (input) INTEGER
           The leading dimension of A.  LDA &gt;= max(1,N).

   B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
           On entry, the second of the pair of matrices whose
           generalized eigenvalues and (optionally) generalized
           eigenvectors are to be computed.
           On exit, the contents will have been destroyed.  (For a
           description of the contents of B on exit, see \"Further
           Details\", below.)

   LDB     (input) INTEGER
           The leading dimension of B.  LDB &gt;= max(1,N).

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
           If JOBVL = 'V', the left generalized eigenvectors.  (See
           \"Purpose\", above.)  Real eigenvectors take one column,
           complex take two columns, the first for the real part and
           the second for the imaginary part.  Complex eigenvectors
           correspond to an eigenvalue with positive imaginary part.
           Each eigenvector will be scaled so the largest component
           will have abs(real part) + abs(imag. part) = 1, *except*
           that for eigenvalues with alpha=beta=0, a zero vector will
           be returned as the corresponding eigenvector.
           Not referenced if JOBVL = 'N'.

   LDVL    (input) INTEGER
           The leading dimension of the matrix VL. LDVL &gt;= 1, and
           if JOBVL = 'V', LDVL &gt;= N.

   VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
           If JOBVR = 'V', the right generalized eigenvectors.  (See
           \"Purpose\", above.)  Real eigenvectors take one column,
           complex take two columns, the first for the real part and
           the second for the imaginary part.  Complex eigenvectors
           correspond to an eigenvalue with positive imaginary part.
           Each eigenvector will be scaled so the largest component
           will have abs(real part) + abs(imag. part) = 1, *except*
           that for eigenvalues with alpha=beta=0, a zero vector will
           be returned as the corresponding eigenvector.
           Not referenced if JOBVR = 'N'.

   LDVR    (input) INTEGER
           The leading dimension of the matrix VR. LDVR &gt;= 1, and
           if JOBVR = 'V', LDVR &gt;= N.

   WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

   LWORK   (input) INTEGER
           The dimension of the array WORK.  LWORK &gt;= max(1,8*N).
           For good performance, LWORK must generally be larger.
           To compute the optimal value of LWORK, call ILAENV to get
           blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:
           NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR;
           The optimal LWORK is:
               2*N + MAX( 6*N, N*(NB+1) ).

           If LWORK = -1, then a workspace query is assumed; the routine
           only calculates the optimal size of the WORK array, returns
           this value as the first entry of the WORK array, and no error
           message related to LWORK is issued by XERBLA.

   INFO    (output) INTEGER
           = 0:  successful exit
           &lt; 0:  if INFO = -i, the i-th argument had an illegal value.
           = 1,...,N:
                 The QZ iteration failed.  No eigenvectors have been
                 calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
                 should be correct for j=INFO+1,...,N.
           &gt; N:  errors that usually indicate LAPACK problems:
                 =N+1: error return from DGGBAL
                 =N+2: error return from DGEQRF
                 =N+3: error return from DORMQR
                 =N+4: error return from DORGQR
                 =N+5: error return from DGGHRD
                 =N+6: error return from DHGEQZ (other than failed
                                                 iteration)
                 =N+7: error return from DTGEVC
                 =N+8: error return from DGGBAK (computing VL)
                 =N+9: error return from DGGBAK (computing VR)
                 =N+10: error return from DLASCL (various calls)

   Further Details
   ===============

   Balancing
   ---------

   This driver calls DGGBAL to both permute and scale rows and columns
   of A and B.  The permutations PL and PR are chosen so that PL*A*PR
   and PL*B*R will be upper triangular except for the diagonal blocks
   A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
   possible.  The diagonal scaling matrices DL and DR are chosen so
   that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to
   one (except for the elements that start out zero.)

   After the eigenvalues and eigenvectors of the balanced matrices
   have been computed, DGGBAK transforms the eigenvectors back to what
   they would have been (in perfect arithmetic) if they had not been
   balanced.

   Contents of A and B on Exit
   -------- -- - --- - -- ----

   If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
   both), then on exit the arrays A and B will contain the real Schur
   form[*] of the \"balanced\" versions of A and B.  If no eigenvectors
   are computed, then only the diagonal blocks will be correct.

   [*] See DHGEQZ, DGEGS, or read the book \"Matrix Computations\",
       by Golub &amp; van Loan, pub. by Johns Hopkins U. Press.

   =====================================================================
</pre>
</html>"));
      end dgegv;

      function dgeqp3 "Obsolete: Computes a QR factorization with column pivoting"

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

        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dgeqp3 instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dgeqp3\">Modelica.Math.Matrices.LAPACK.dgeqp3</a>
instead.
</p>
<p>
Note: regarding the function from the Modelica Standard Library 
<p>
<ul>
  <li>
    outputs <code>jpvt</code> and <code>tau</code> are interconverted,
  </li>
  <li>
    output <code>work</code> is missing.
  </li>
</ul>
<p>
This means a function call like
</p>
<blockquote><pre>
(Aout, jpvt, tau, info, work) :=
  ObsoleteLinearSystems2.Math.LAPACK.dgeqp3(A, lwork1);
</pre></blockquote>
<p>
could be replaced with
</p>
<blockquote><pre>
(Aout, tau, jpvt, info) := 
  Modelica.Math.Matrices.LAPACK.dgeqp3(A, lwork);
</pre></blockquote>
</html>"));
      end dgeqp3;

      function dgeqrf "Obsolete: Computes a QR factorization without pivoting"

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

        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dgeqrf instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dgeqrf\">Modelica.Math.Matrices.LAPACK.dgeqrf</a>
instead.
</p>
<p>
Note: input <code>lwork1</code> is missing in the function from the Modelica
Standard Library and must be removed or resolved in another way.
</p>
</html>"));
      end dgeqrf;

      function dgetrs
        "Obsolete: Solves a system of linear equations with the LU decomposition from dgetrf(..)"

        input Real LU[:,size(LU, 1)] "LU factorization of dgetrf of a square matrix";
        input Integer pivots[size(LU, 1)] "Pivot vector of dgetrf";
        input Real B[size(LU, 1),:] "Right hand side matrix B";
        output Real X[size(B, 1),size(B, 2)]=B "Solution matrix X";

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
          info) annotation (Library="lapack");
        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dgetrs instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dgetrs\">Modelica.Math.Matrices.LAPACK.dgetrs</a>
instead.
</p>
<p>
Note: additional output <code>info</code> exists in the function from the Modelica
Standard Library. It must be added or resolved in another way.
</p>
</html>"));
      end dgetrs;

      function dhseqr
        "Obsolete: Compute eingenvalues of a matrix A using lapack routine DHSEQR for Hessenberg form matrix"
        input Real H[:,size(H, 1)];
        input Integer lwork=max(1, size(H, 1));
        input Boolean eigenValuesOnly=true;
        input String compz="N";
        input Real Z[:,:]=H;
        output Real alphaReal[size(H, 1)]
          "Real part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
        output Real alphaImag[size(H, 1)]
          "Imaginary part of alpha (eigenvalue=(alphaReal+i*alphaImag))";
        output Integer info;
        output Real Ho[:,:]=H;
        output Real Zo[:,:]=Z;
        output Real work[max({lwork, size(H, 1),1})];

      protected
        Integer n=size(H, 1);
        String job=if eigenValuesOnly then "E" else "S";
        Integer ilo=1;
        Integer ihi=n;
        Integer ldh=max(n, 1);
        Integer lw=if lwork == -1 then -1 else max(lwork, size(H, 1));

      external "Fortran 77" dhseqr(
          job,
          compz,
          n,
          ilo,
          ihi,
          Ho,
          ldh,
          alphaReal,
          alphaImag,
          Zo,
          ldh,
          work,
          lw,
          info) annotation(Library = {"lapack"});

        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.LAPACK.dhseqr instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.LAPACK.dhseqr\">Modelica.Math.Matrices.LAPACK.dhseqr</a>
instead.
</p>
<p>
Note: input <code>lwork</code> is missing in the function from the Modelica
Standard Library and must be removed or resolved in another way.
</p>
</html>"));
      end dhseqr;
      annotation (Documentation(info="<html>
<p>
This package contains functions to call routines from software library
<a href=\"https://www.netlib.org/lapack/\">LAPACK </a>
(Linear Algebra PACKage) aimed for numerical linear algebra. The library is
provided by <a href=\"https://www.netlib.org/\">Netlib Repository</a>.
</p>
</html>"));
    end LAPACK;

    package Matrices "Package of functions operating on matrices"
      extends Modelica.Icons.Package;

      function leastSquares "Obsolete: Solve overdetermined or underdetermined real system of linear equations A*x=b in a least squares sense (A may be rank deficient)"
        extends Modelica.Icons.Function;
        input Real A[:,:] "Matrix A";
        input Real b[size(A, 1)] "Vector b";
        output Real x[size(A, 2)]
          "Vector x such that min|A*x-b|^2 if size(A,1) >= size(A,2) or min|x|^2 and A*x=b, if size(A,1) < size(A,2)";

      protected
        Integer rank;
      algorithm
        (x,rank) :=Modelica.Math.Matrices.leastSquares(A,b);
        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Matrices.leastSquares instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Matrices.leastSquares\">Modelica.Math.Matrices.leastSquares</a>
instead.
</p>
</html>"));
      end leastSquares;
    end Matrices;

    package Internal "Package of internal functions operating on matrices (for advanced users only)"
      extends Modelica.Icons.InternalPackage;


      function dhseqr_workdim "Calculate the optimal size of the WORK array in dhseqr"
        import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input Real H[:,:];

        output Integer lwork;
        output Integer info;

      protected
        Real work[:];

      algorithm
        if min(size(H, 1), size(H, 2)) > 0 then
          (,,info,,,work) := ObsoleteLinearSystems2.Math.LAPACK.dhseqr(H, -1);
          lwork := integer(work[1]);
        else
          lwork := 1;
        end if;

        annotation (
          obsolete = "Deprecated function - see ObsoleteLinearSystems2.Math.LAPACK.dhseqr for alternatives");
      end dhseqr_workdim;

      function dgeqrf_workdim "Calculate the optimal size of the WORK array in dgeqrf"
        import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input Real A[:,:];

        output Integer lwork;
        output Integer info;

      protected
        Real work[max(1,size(A,1))];

      algorithm
        (,,info,work) := ObsoleteLinearSystems2.Math.LAPACK.dgeqrf(A, -1);
        lwork := integer(work[1]);

        annotation (
          obsolete = "Deprecated function - see ObsoleteLinearSystems2.Math.LAPACK.dgeqrf for alternatives");
      end dgeqrf_workdim;
      annotation (Documentation(info="<html>
<p>
Generally, the functions in this package should not be used by the user.
</p>
<p>
This package contains functions which cannot be used in an arbitrary
way and require particular knowledge.
Therefore, only advanced users should deal with contained classes.
</p>
</html>"));
    end Internal;
  end Math;

  annotation (
    preferredView="info",
    version="3.0.0",
    versionDate="2024-06-21",
    dateModified = "2024-03-26 10:00:00Z",
    revisionId="$Format:%h %ci$",
    uses(
      Modelica(version="4.0.0"),
      Complex(version="4.0.0")),
    Documentation(info="<html>
<p>
This package contains functions and blocks from the Modelica_LinearSystems2 Library
version 2.4.X that are no longer available in the version 3.0.0.
The conversion script for version 2.4.X changes references in existing
user models automatically to the functions and blocks of package
ObsoleteLinearSystems2. The user should <strong>manually</strong> replace all
references to ObsoleteLinearSystems2 in his/her models to the functions
and blocks that are recommended in the documentation of the respective model.
</p>

<p>
In most cases, this means that a&nbsp;model with the name
\"ObsoleteLinearSystems2.XY\" should be renamed to \"Modelica_LinearSystems2.YZ\"
(version 3.0.0) and manually adaptated afterwards.
This usually requires some changes at the place where
the class is used (besides the renaming of the underlying class).
</p>

<p>
The models in ObsoleteLinearSystems2 are either not according to the
Modelica Language version 3.6 and higher, or the model was changed to get
a&nbsp;better design.
In all cases, an automatic conversion to the new implementation
was not feasible, since too complicated.
See also 
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.ReleaseNotes.Version_3_0_0\">Modelica_LinearSystems2.UsersGuide.ReleaseNotes.Version_3_0_0</a>.
</p>

<p>
In order to easily detect obsolete models and blocks, all of them are specially
marked in the icon layer with a red box. Additionally, an annotation &quot;obsolete&quot;
is provided.
</p>

<p>
<strong>Copyright</strong> &copy; 2024, DLR Institute of System Dynamics and Control
</p>

<p>
<em>
This Modelica package is <u>free</u> software and
the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the
3-Clause BSD license, see the license conditions (including the
disclaimer of warranty) in the
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.The3clauseBSDLicense\">User's Guide</a>.
</em>
</p>
</html>"));
end ObsoleteLinearSystems2;
