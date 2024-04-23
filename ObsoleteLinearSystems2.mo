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

    package Vectors "Package of functions operating on vectors"
      extends Modelica.Icons.Package;

      function printVector "Print vector"
        import Modelica.Utilities.Strings;

        input Real v[:] "Vector of real numbers to be printed";
        input Integer significantDigits=6
          "Number of significant digits that are shown";
        input String name="v" "Independent variable name used for printing";
        output String s="" "String containing v";
      protected
        String blanks=Strings.repeat(significantDigits);
        String space=Strings.repeat(8);
        String space2=Strings.repeat(3);
        Integer r=size(v, 1);

      algorithm
        if r == 0 then
          s := name + " = []";
        else
          s := "\n" + name + " = \n";
          for i in 1:r loop
            s := s + space;

            if v[i] >= 0 then
              s := s + " ";
            end if;
            s := s + String(v[i], significantDigits=significantDigits) +
              Strings.repeat(significantDigits + 8 - Strings.length(String(abs(v[i]))));

            s := s + "\n";
          end for;

        end if;
        annotation (
          obsolete = "Deprecated function - use Modelica.Math.Vectors.toString instead",
          Documentation(info="<html>
<p>
This function is obsolete. Use
<a href=\"modelica://Modelica.Math.Vectors.toString\">Modelica.Math.Vectors.toString</a>
instead.
</p>
<p>
Note: the inputs two and three (<code>significantDigits</code> and <code>name</code>) are
interchanged in Modelica.Utilities.Strings.isEqual. Consequently, a call like
</p>
<blockquote><pre>
ObsoleteLinearSystems2.Math.Vectors.printVector({3,33,7}, 2, \"vec\");
</pre></blockquote>
<p>
shall be replaced with
</p>
<blockquote><pre>
Modelica.Math.Vectors.toString({3,33,7}, \"vec\", 2);
</pre></blockquote>
</html>"));
      end printVector;
    end Vectors;

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

  package Controller "Package of continuous and discrete input/output blocks"
    extends Modelica.Icons.Package;

    package Examples
      extends Modelica.Icons.ExamplesPackage;

      package Components
        extends Modelica.Icons.UtilitiesPackage;
        model DoublePendulum "crane trolley system"
          extends Modelica.Icons.ObsoleteModel;

          parameter Modelica.Units.SI.Mass m_trolley = 5;
          parameter Modelica.Units.SI.Mass m_load = 20;
          parameter Modelica.Units.SI.Length length = 2;
          parameter Modelica.Units.SI.Angle phi1_start = -80.0/180*pi;
          parameter Modelica.Units.SI.Angle phi2_start = 10;
          parameter Modelica.Units.SI.AngularVelocity w1_start = 0.0;
          parameter Modelica.Units.SI.AngularVelocity w2_start = 0.0;

          constant Real pi = Modelica.Constants.pi;

          inner Modelica.Mechanics.MultiBody.World world(animateWorld=false,
              animateGravity=false)
            annotation (Placement(transformation(extent={{-140,-80},{-120,-60}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
            annotation (Placement(transformation(extent={{-96,0},{-76,20}})));
          Modelica.Mechanics.Translational.Components.Damper damper1(d=0)
            annotation (Placement(transformation(extent={{-96,14},{-76,34}})));
          Modelica.Mechanics.MultiBody.Joints.Revolute rev(n={0,0,1},useAxisFlange=true,
            phi(fixed=true, start=phi1_start),
            w(fixed=true, start=w1_start))
            annotation (Placement(transformation(extent={{-30,0},{-10,20}}, rotation=0)));
          Modelica.Mechanics.Rotational.Components.Damper damper(d=0)
            annotation (Placement(transformation(extent={{-22,40},{-2,60}},rotation=0)));
          Modelica.Mechanics.MultiBody.Parts.Body body(
            m=m_load,
            r_CM={0,0,0},
            specularCoefficient=4*world.defaultSpecularCoefficient,
            sphereDiameter=1.5*world.defaultBodyDiameter)
            annotation (Placement(transformation(extent={{78,0},{98,20}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
            shapeType="box",
            m=m_trolley,
            sphereDiameter=world.defaultBodyDiameter,
            r={0,0,0},
            r_CM={0,0,0})
            annotation (Placement(transformation(extent={{-58,-2},{-38,18}})));
          Modelica.Mechanics.Translational.Sources.Force force
            annotation (Placement(transformation(extent={{-98,34},{-78,54}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles
            annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeVelocity relativeVelocity
            annotation (Placement(transformation(extent={{-96,-30},{-76,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition
            annotation (Placement(transformation(extent={{-96,-60},{-76,-40}})));
          Modelica.Blocks.Interfaces.RealInput u
            annotation (Placement(transformation(extent={{-190,-20},{-150,20}})));
          Modelica.Blocks.Interfaces.RealOutput s
            annotation (Placement(transformation(extent={{150,90},{170,110}})));
          Modelica.Blocks.Interfaces.RealOutput v
            annotation (Placement(transformation(extent={{150,50},{170,70}})));
         Modelica.Blocks.Interfaces.RealOutput phi
            annotation (Placement(transformation(extent={{150,10},{170,30}})));
          Modelica.Blocks.Interfaces.RealOutput w
            annotation (Placement(transformation(extent={{150,-30},{170,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity
            annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

          Modelica.Blocks.Sources.Constant const(k=0.5*Modelica.Constants.pi)
            annotation (Placement(transformation(extent={{94,-22},{106,-10}})));
          Modelica.Blocks.Math.Add add
            annotation (Placement(transformation(extent={{116,-10},{136,10}})));
          Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
            phi(fixed=true, start=phi2_start),
            w(fixed=true, start=w2_start),
            cylinderDiameter=3*world.defaultJointWidth,
            cylinderColor={0,0,200}) annotation (Placement(transformation(extent={{24,0},{
                    44,20}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles1
            annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity1
            annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
         Modelica.Blocks.Interfaces.RealOutput phi1
            annotation (Placement(transformation(extent={{150,-70},{170,-50}})));
          Modelica.Blocks.Interfaces.RealOutput w1
            annotation (Placement(transformation(extent={{150,-110},{170,-90}})));
          Modelica.Blocks.Math.Add add1
            annotation (Placement(transformation(extent={{88,-50},{108,-30}})));
          Modelica.Blocks.Sources.Constant const1(k=0)
            annotation (Placement(transformation(extent={{66,-62},{78,-50}})));
          Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
            r={length/2,0,0},
            specularCoefficient=0.7,
            color={0,0,0},
            diameter=0.05,
            density=900)
            annotation (Placement(transformation(extent={{-4,0},{16,20}})));
          Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
            r={length/2,0,0},
            specularCoefficient=0.7,
            color={0,0,0},
            diameter=0.05,
            density=900)
            annotation (Placement(transformation(extent={{52,0},{72,20}})));
        equation
          connect(damper.flange_b, rev.axis) annotation (Line(points={{-2,50},{0,50},{0,
                  24},{0,20},{-20,20}}, color={0,0,0}));
          connect(rev.support, damper.flange_a) annotation (Line(points={{-26,20},{-26,
                  26},{-36,26},{-36,50},{-22,50}}, color={0,0,0}));
          connect(bodyShape.frame_b, rev.frame_a) annotation (Line(
              points={{-38,8},{-34,8},{-34,10},{-30,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(prismatic.frame_a, world.frame_b) annotation (Line(
              points={{-96,10},{-110,10},{-110,-70},{-120,-70}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(force.flange, prismatic.axis) annotation (Line(
              points={{-78,44},{-78,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(damper1.flange_a, prismatic.support) annotation (Line(
              points={{-96,24},{-96,16},{-90,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(damper1.flange_b, prismatic.axis) annotation (Line(
              points={{-76,24},{-78,24},{-78,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
              points={{-76,10},{-68,10},{-68,8},{-58,8}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeVelocity.frame_b, prismatic.frame_b) annotation (Line(
              points={{-76,-20},{-76,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeVelocity.frame_a, prismatic.frame_a) annotation (Line(
              points={{-96,-20},{-96,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativePosition.frame_b, relativeVelocity.frame_b) annotation (Line(
              points={{-76,-50},{-76,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativePosition.frame_a, relativeVelocity.frame_a) annotation (Line(
              points={{-96,-50},{-96,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles.frame_b, rev.frame_b) annotation (Line(
              points={{-10,-20},{-10,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles.frame_a, rev.frame_a) annotation (Line(
              points={{-30,-20},{-30,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(u, force.f) annotation (Line(
              points={{-170,0},{-136,0},{-136,44},{-100,44}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngularVelocity.frame_a, relativeAngles.frame_a) annotation (
              Line(
              points={{-30,-50},{-30,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity.frame_b, relativeAngles.frame_b) annotation (
              Line(
              points={{-10,-50},{-10,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity.w_rel[3], w) annotation (Line(
              points={{-20,-60.6667},{-20,-66},{120,-66},{120,-20},{160,-20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeVelocity.v_rel[1], v) annotation (Line(
              points={{-86,-31.3333},{-104,-31.3333},{-104,-32},{-118,-32},{-118,62},{42,62},{42,60},{160,60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativePosition.r_rel[1], s) annotation (Line(
              points={{-86,-61.3333},{-104,-61.3333},{-104,-58},{-120,-58},{-120,100},{160,100}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(phi, phi) annotation (Line(
              points={{160,20},{160,20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add.y, phi) annotation (Line(
              points={{137,0},{148,0},{148,20},{160,20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(const.y, add.u2) annotation (Line(
              points={{106.6,-16},{110,-16},{110,-6},{114,-6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add.u1, relativeAngles.angles[3]) annotation (Line(
              points={{114,6},{108,6},{108,-4},{58,-4},{58,-36},{-20,-36},{-20,-30.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngles1.frame_a, revolute2.frame_a) annotation (Line(
              points={{24,-20},{24,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles1.frame_b, revolute2.frame_b) annotation (Line(
              points={{44,-20},{44,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles1.frame_a, relativeAngularVelocity1.frame_a)
            annotation (Line(
              points={{24,-20},{24,-50}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity1.frame_b, relativeAngles1.frame_b)
            annotation (Line(
              points={{44,-50},{44,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(const1.y, add1.u2)
            annotation (Line(
              points={{78.6,-56},{82,-56},{82,-46},{86,-46}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add1.u1, relativeAngles1.angles[3]) annotation (Line(
              points={{86,-34},{60,-34},{60,-30.6667},{34,-30.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add1.y, phi1) annotation (Line(
              points={{109,-40},{136,-40},{136,-60},{160,-60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngularVelocity1.w_rel[3], w1) annotation (Line(
              points={{34,-60.6667},{36,-60.6667},{36,-100},{160,-100}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
              points={{72,10},{78,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
              points={{52,10},{44,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
              points={{16,10},{24,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder.frame_a, rev.frame_b) annotation (Line(
              points={{-4,10},{-10,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (
            obsolete = "Deprecated model - use Modelica_LinearSystems2.Utilities.Plants.DoublePendulum instead",
            Diagram(coordinateSystem(
                preserveAspectRatio=true,
                extent={{-150,-100},{150,100}},
                grid={2,2}), graphics),
            Documentation(info="<html>
<p>
This plant model is obsolete. Use
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\">Utilities.Plants.DoublePendulum</a>
instead.
</p>
</html>"),  Icon(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{150,
                    100}}), graphics={
                Rectangle(
                  extent={{-150,122},{150,-120}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-82,22},{82,18}},
                  lineColor={0,0,255},
                  fillPattern=FillPattern.Forward),
                Rectangle(extent={{-44,54},{0,28}}, lineColor={0,0,0}),
                Ellipse(
                  extent={{-40,34},{-28,22}},
                  lineColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,255,255},
                  lineThickness=0.5),
                Ellipse(
                  extent={{-16,34},{-4,22}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineThickness=0.5),
                Line(
                  points={{-18,-16},{10,-62}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Ellipse(
                  extent={{4,-56},{20,-72}},
                  lineColor={0,0,0},
                  fillColor={0,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-25,44},{-19,38}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),
                Line(
                  points={{28,46},{4,46}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{34,40},{10,40}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{-22,40},{-18,-16}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Ellipse(
                  extent={{-20,-15},{-14,-21}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),Rectangle(
                  extent={{-152,124},{152,-122}},
                  lineColor={255,0,0},
                  pattern=LinePattern.Dash,
                  lineThickness=0.5)}));
        end DoublePendulum;

        model InverseDoublePendulum "Inverse double pendulum"
          extends Modelica.Icons.ObsoleteModel;

          parameter Modelica.Units.SI.Mass m_trolley = 1 "Mass of trolley";
          parameter Modelica.Units.SI.Mass m_load = 1 "Mass of load on 2nd arm";
          parameter Modelica.Units.SI.Length length = 1
            "Total length of double pendulum (i.e. length of each arm = length/2)";
          parameter Modelica.Units.SI.Angle phi1_start = 90.0/180*pi
            "Initial rotation angle of 1st arm relative to trolley";
          parameter Modelica.Units.SI.Angle phi2_start = 0
            "Initial rotation angle of 2nd arm relative to 1st arm";
          parameter Modelica.Units.SI.AngularVelocity w1_start = 0.0
            "Initial angular velocity of 1st arm relative to trolley";
          parameter Modelica.Units.SI.AngularVelocity w2_start = 0.0
            "Initial angular velocity of 2nd arm relative to 1st arm";

          parameter Modelica.Units.SI.Position s_start = 0.0
            "Initial position of trolley relative to world";
          parameter Modelica.Units.SI.Velocity v_start = 0.0
            "Initial velocity of trolley relative to world";

          parameter Boolean cartDisturbance=false
            "True, if cart disturbance should be enabled";
          parameter Boolean bodyDisturbance=false
            "True, if body disturbance should be enabled";

          constant Real pi=Modelica.Constants.pi;

          inner Modelica.Mechanics.MultiBody.World world(
            gravityType=Modelica.Mechanics.MultiBody.Types.GravityTypes.UniformGravity,
            animateWorld=false,
            animateGravity=false) annotation (Placement(transformation(extent={{-140,-80},
                    {-120,-60}}, rotation=0)));

          Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(
            useAxisFlange=true,
            s(start=s_start, fixed=true),
            animation=false,
            v(start=v_start, fixed=true))
            annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
          Modelica.Mechanics.Translational.Components.Damper damper1(d=0)
            annotation (Placement(transformation(extent={{-100,18},{-80,38}})));
          Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(
            n={0,0,1},
            useAxisFlange=true,
            phi(fixed=true, start=phi1_start),
            w(fixed=true, start=w1_start),
            cylinderColor=bodyShape.color,
            cylinderDiameter=2*bodyCylinder.diameter) annotation (Placement(
                transformation(extent={{-30,0},{-10,20}}, rotation=0)));
          Modelica.Mechanics.Rotational.Components.Damper damper(d=0) annotation (
              Placement(transformation(extent={{-30,22},{-10,42}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Parts.Body body(
            m=m_load,
            r_CM={0,0,0},
            specularCoefficient=4*world.defaultSpecularCoefficient,
            sphereDiameter=1.5*world.defaultBodyDiameter,
            sphereColor=bodyCylinder1.color) annotation (Placement(transformation(
                  extent={{80,0},{100,20}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
            m=m_trolley,
            sphereDiameter=world.defaultBodyDiameter,
            animateSphere=false,
            shapeType="box",
            lengthDirection={0,-1,0},
            widthDirection={1,0,0},
            length=0.1,
            r_shape=0.3*revolute1.cylinderDiameter*{0,-1,0},
            width=0.3,
            height=0.5*bodyShape.width,
            color={0,0,0},
            r={0,0,0},
            r_CM={0,0,0})
            annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
          Modelica.Mechanics.Translational.Sources.Force force
            annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles
            annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeVelocity relativeVelocity
            annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition
            annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity
            annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

          Modelica.Blocks.Sources.Constant const(k=-0.5*Modelica.Constants.pi)
            annotation (Placement(transformation(extent={{98,34},{110,46}})));
          Modelica.Blocks.Math.Add add
            annotation (Placement(transformation(extent={{120,30},{140,10}})));
          Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
            phi(fixed=true, start=phi2_start),
            w(fixed=true, start=w2_start),
            specularCoefficient=0.7,
            cylinderDiameter=2*bodyCylinder.diameter,
            cylinderColor=bodyCylinder.color)
            annotation (Placement(transformation(extent={{24,0},{44,20}}, rotation=0)));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles1
            annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
          Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity1
            annotation (Placement(transformation(extent={{24,-64},{44,-44}})));
          Modelica.Blocks.Sources.Constant const1(k=0)
            annotation (Placement(transformation(extent={{98,-82},{110,-70}})));
          Modelica.Blocks.Math.Add add1
            annotation (Placement(transformation(extent={{120,-70},{140,-50}})));
          Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                  extent={{-190,-20},{-150,20}}), iconTransformation(extent={{-140,-20},
                    {-100,20}})));
          Modelica.Blocks.Interfaces.RealOutput s annotation (Placement(transformation(
                  extent={{150,90},{170,110}}), iconTransformation(extent={{100,90},{
                    120,110}})));
          Modelica.Blocks.Interfaces.RealOutput v annotation (Placement(transformation(
                  extent={{150,50},{170,70}}), iconTransformation(extent={{100,50},{120,
                    70}})));
          Modelica.Blocks.Interfaces.RealOutput phi annotation (Placement(
                transformation(extent={{150,10},{170,30}}), iconTransformation(extent={
                    {100,10},{120,30}})));
          Modelica.Blocks.Interfaces.RealOutput w annotation (Placement(transformation(
                  extent={{150,-30},{170,-10}}), iconTransformation(extent={{100,-30},{
                    120,-10}})));
          Modelica.Blocks.Interfaces.RealOutput phi1 annotation (Placement(
                transformation(extent={{150,-70},{170,-50}}), iconTransformation(extent=
                   {{100,-70},{120,-50}})));
          Modelica.Blocks.Interfaces.RealOutput w1 annotation (Placement(transformation(
                  extent={{150,-110},{170,-90}}), iconTransformation(extent={{100,-110},
                    {120,-90}})));
          Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
            r={length/2,0,0},
            specularCoefficient=0.7,
            diameter=0.03,
            density=1000,
            color={0,128,0})
            annotation (Placement(transformation(extent={{-2,0},{18,20}})));
          Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
            r={length/2,0,0},
            specularCoefficient=0.7,
            diameter=0.03,
            density=1000,
            color={0,0,255})
            annotation (Placement(transformation(extent={{52,0},{72,20}})));
          Modelica.Blocks.Interfaces.RealInput dist if cartDisturbance annotation (
              Placement(transformation(
                extent={{-20,-20},{20,20}},
                rotation=-90,
                origin={-80,120}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=-90,
                origin={-60,120})));
          Modelica.Mechanics.Translational.Sources.Force distrubanceForceCart if
            cartDisturbance annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=-90,
                origin={-80,80})));
          Modelica.Blocks.Interfaces.RealInput dist2 if bodyDisturbance annotation (
              Placement(transformation(
                extent={{-20,-20},{20,20}},
                rotation=-90,
                origin={80,120}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=-90,
                origin={60,120})));
          Modelica.Mechanics.MultiBody.Forces.Torque torque if bodyDisturbance
            annotation (Placement(transformation(extent={{40,60},{60,80}})));
          Modelica.Blocks.Sources.Constant const2[2](k={0,0}) if bodyDisturbance
            annotation (Placement(transformation(extent={{0,70},{20,90}})));
          Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape1(
            shapeType="cylinder",
            lengthDirection=revolute2.n,
            widthDirection={0,0,1},
            length=1.2*revolute2.cylinderLength,
            width=bodyCylinder1.diameter,
            height=bodyCylinder1.diameter,
            color=bodyCylinder1.color,
            r_shape=-0.5*fixedShape1.length*fixedShape1.lengthDirection)
            annotation (Placement(transformation(extent={{52,20},{72,40}})));
          Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
            shapeType="cylinder",
            lengthDirection=revolute1.n,
            widthDirection={0,0,1},
            length=1.2*revolute1.cylinderLength,
            width=bodyCylinder.diameter,
            height=bodyCylinder.diameter,
            color=bodyCylinder.color,
            r_shape=-0.5*fixedShape.length*fixedShape.lengthDirection)
            annotation (Placement(transformation(extent={{0,20},{20,40}})));
          Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape2(
            lengthDirection={0,-1,0},
            widthDirection={1,0,0},
            height=0.5*bodyShape.height,
            color={100,100,100},
            length=0.5*bodyShape.length,
            width=50*bodyShape.width,
            r_shape=bodyShape.r_shape + 0.5*(bodyShape.length - fixedShape2.length)*{
                200,-1,0})
            annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
        equation
          connect(damper.flange_b, revolute1.axis) annotation (Line(points={{-10,32},{-8,
                  32},{-8,22},{-8,20},{-20,20}}, color={0,0,0}));
          connect(revolute1.support, damper.flange_a) annotation (Line(points={{-26,20},
                  {-26,20},{-34,20},{-34,32},{-30,32}}, color={0,0,0}));
          connect(bodyShape.frame_b, revolute1.frame_a) annotation (Line(
              points={{-38,10},{-30,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(prismatic.frame_a, world.frame_b) annotation (Line(
              points={{-100,10},{-110,10},{-110,-70},{-120,-70}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(force.flange, prismatic.axis) annotation (Line(
              points={{-80,50},{-70,50},{-70,16},{-82,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(damper1.flange_a, prismatic.support) annotation (Line(
              points={{-100,28},{-100,16},{-94,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(damper1.flange_b, prismatic.axis) annotation (Line(
              points={{-80,28},{-80,16},{-82,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
              points={{-80,10},{-58,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeVelocity.frame_b, prismatic.frame_b) annotation (Line(
              points={{-80,-20},{-80,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeVelocity.frame_a, prismatic.frame_a) annotation (Line(
              points={{-100,-20},{-100,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativePosition.frame_b, relativeVelocity.frame_b) annotation (Line(
              points={{-80,-50},{-80,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativePosition.frame_a, relativeVelocity.frame_a) annotation (Line(
              points={{-100,-50},{-100,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles.frame_b, revolute1.frame_b) annotation (Line(
              points={{-10,-20},{-10,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles.frame_a, revolute1.frame_a) annotation (Line(
              points={{-30,-20},{-30,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(u, force.f) annotation (Line(
              points={{-170,0},{-136,0},{-136,50},{-102,50}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngularVelocity.frame_a, relativeAngles.frame_a) annotation (
              Line(
              points={{-30,-50},{-30,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity.frame_b, relativeAngles.frame_b) annotation (
              Line(
              points={{-10,-50},{-10,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity.w_rel[3], w) annotation (Line(
              points={{-20,-60.6667},{-20,-72},{80,-72},{80,-30},{140,-30},{140,-20},{160,-20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeVelocity.v_rel[1], v) annotation (Line(
              points={{-90,-31.3333},{-108,-31.3333},{-108,-34},{-122,-34},{-122,60},{160,60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativePosition.r_rel[1], s) annotation (Line(
              points={{-90,-61.3333},{-108,-61.3333},{-108,-58},{-124,-58},{-124,100},{160,100}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add.y, phi) annotation (Line(
              points={{141,20},{160,20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(const.y, add.u2) annotation (Line(
              points={{110.6,40},{114,40},{114,26},{118,26}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add.u1, relativeAngles.angles[3]) annotation (Line(
              points={{118,14},{114,14},{114,-20},{60,-20},{60,-36},{-20,-36},{-20,-30.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngles1.frame_a, revolute2.frame_a) annotation (Line(
              points={{24,-20},{24,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles1.frame_b, revolute2.frame_b) annotation (Line(
              points={{44,-20},{44,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngles1.frame_a, relativeAngularVelocity1.frame_a)
            annotation (Line(
              points={{24,-20},{24,-54}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(relativeAngularVelocity1.frame_b, relativeAngles1.frame_b)
            annotation (Line(
              points={{44,-54},{44,-20}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(const1.y, add1.u2) annotation (Line(
              points={{110.6,-76},{114,-76},{114,-66},{118,-66}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add1.u1, relativeAngles1.angles[3]) annotation (Line(
              points={{118,-54},{70,-54},{70,-34},{34,-34},{34,-30.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add1.y, phi1) annotation (Line(
              points={{141,-60},{160,-60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(relativeAngularVelocity1.w_rel[3], w1) annotation (Line(
              points={{34,-64.6667},{34,-100},{160,-100}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
              points={{18,10},{24,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder.frame_a, revolute1.frame_b) annotation (Line(
              points={{-2,10},{-10,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
              points={{72,10},{80,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
              points={{52,10},{44,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(distrubanceForceCart.f, dist) annotation (Line(
              points={{-80,92},{-80,120}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(distrubanceForceCart.flange, prismatic.axis) annotation (Line(
              points={{-80,70},{-70,70},{-70,16},{-82,16}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(dist2, torque.torque[3]) annotation (Line(
              points={{80,120},{80,90},{44,90},{44,82.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(const2.y, torque.torque[1:2]) annotation (Line(
              points={{21,80},{34,80},{34,86},{44,86},{44,82}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(torque.frame_a, revolute2.frame_a) annotation (Line(
              points={{40,70},{40,70},{24,70},{24,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(torque.frame_b, bodyCylinder1.frame_b) annotation (Line(
              points={{60,70},{74,70},{74,10},{72,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(fixedShape1.frame_a, revolute2.frame_b) annotation (Line(
              points={{52,30},{48,30},{48,10},{44,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(fixedShape.frame_a, revolute1.frame_b) annotation (Line(
              points={{0,30},{-6,30},{-6,10},{-10,10}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          connect(fixedShape2.frame_a, world.frame_b) annotation (Line(
              points={{-100,-90},{-110,-90},{-110,-70},{-120,-70}},
              color={95,95,95},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (
            obsolete = "Deprecated model - use Modelica_LinearSystems2.Utilities.Plants.DoublePendulumInverse instead",
            Diagram(coordinateSystem(
                preserveAspectRatio=true,
                extent={{-150,-100},{150,100}},
                grid={2,2}), graphics),
            Documentation(info="<html>
<p>
This plant model is obsolete. Use
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plants.DoublePendulumInverse\">Utilities.Plants.DoublePendulumInverse</a>
instead.
</p>
</html>"),  Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                    100}}), graphics={
                Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-82,-74},{82,-78}},
                  lineColor={0,0,255},
                  fillPattern=FillPattern.Forward),
                Rectangle(extent={{-44,-42},{0,-68}}, lineColor={0,0,0}),
                Ellipse(
                  extent={{-40,-62},{-28,-74}},
                  lineColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,255,255},
                  lineThickness=0.5),
                Ellipse(
                  extent={{-16,-62},{-4,-74}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineThickness=0.5),
                Line(
                  points={{-16,-6},{-22,-56}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  thickness=0.5),
                Ellipse(
                  extent={{-32,54},{-20,42}},
                  lineColor={0,0,0},
                  fillColor={0,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-25,-52},{-19,-58}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),
                Line(
                  points={{34,-56},{10,-56}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{28,-64},{4,-64}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{-34,54},{-34,54},{-36,52},{-36,46},{-34,44}},
                  color={0,0,255},
                  smooth=Smooth.None),
                Line(
                  points={{-38,56},{-38,56},{-40,54},{-40,48},{-38,46}},
                  color={0,0,255},
                  smooth=Smooth.None),
                Line(
                  points={{-26,42},{-16,-6}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  thickness=0.5),
                Ellipse(
                  extent={{-20,-2},{-14,-8}},
                  lineColor={0,0,0},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{150,
                    100}})));
        end InverseDoublePendulum;
      end Components;
    end Examples;
  end Controller;
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
