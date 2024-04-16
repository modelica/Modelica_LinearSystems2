within ;
package ObsoleteLinearSystems2
  "Library that contains components from Modelica_LinearSystems2 Library 2.4.X that have been removed from version 3.0.0"

  package Math "Package of additional functions for Modelica.Math"
    extends Modelica.Icons.Package;

    operator record Complex "Record defining a Complex number"

      encapsulated function j "Obsolete imaginary unit function"
        import .Modelica;
        import .Complex;

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

        annotation (Documentation(info="<html>
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

        annotation (Documentation(info="<html>
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

        annotation (Documentation(info="<html>
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
        annotation (Documentation(info="<html>
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

        annotation (Documentation(info="<html>
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
