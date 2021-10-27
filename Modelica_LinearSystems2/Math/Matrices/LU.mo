within Modelica_LinearSystems2.Math.Matrices;
function LU "LU decomposition of square or rectangular matrix"
  extends Modelica.Icons.Function;
  input Real A[:,:] "Square or rectangular matrix";
  output Real LU[size(A, 1),size(A, 2)]=A
    "L,U factors (used with LU_solve(..))";
  output Integer pivots[min(size(A, 1), size(A, 2))]
    "Pivot indices (used with LU_solve(..))";
  output Integer info "Information";

algorithm
  if min(size(A)) > 0 then
    (LU,pivots,info) := Modelica.Math.Matrices.LAPACK.dgetrf(A);
  else
    LU := fill(
      0,
      size(A, 1),
      size(A, 2));
    pivots := fill(0, min(size(A, 1), size(A, 2)));
    info := 0;
  end if;

  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.LU instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(LU, pivots)       = Matrices.<strong>LU</strong>(A);
(LU, pivots, info) = Matrices.<strong>LU</strong>(A);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
LU decomposition of a &quot;Real[m,n]&quot; matrix A, i.e.,
</p>
<blockquote>
<p>
<strong>P</strong>*<strong>L</strong>*<strong>U</strong> = <strong>A</strong>
</p>
</blockquote>
<p>
where <strong>P</strong> is a permutation matrix (implicitely
defined by vector <code>pivots</code>),
<strong>L</strong> is a lower triangular matrix with unit
diagonal elements (lower trapezoidal if m &gt; n), and
<strong>U</strong> is an upper triangular matrix (upper trapezoidal if m &lt; n).
Matrices <strong>L</strong> and <strong>U</strong> are stored in the returned
matrix <code>LU</code> (the diagonal of <strong>L</strong> is not stored).
With the companion function
<a href=\"modelica://Modelica.Math.Matrices.LU_solve\">Modelica.Math.Matrices.LU_solve</a>,
this decomposition can be used to solve
linear systems (<strong>P</strong>*<strong>L</strong>*<strong>U</strong>)*<strong>x</strong> = <strong>b</strong> with different right
hand side vectors <strong>b</strong>. If a linear system of equations with
just one right hand side vector <strong>b</strong> shall be solved, it is
more convenient to just use the function
<a href=\"modelica://Modelica.Math.Matrices.solve\">Modelica.Math.Matrices.solve</a>.
</p>
<p>
The optional third (Integer) output argument has the following meaning:
</p>
<blockquote>
<table border=0 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\">info = 0:</td>
      <td valign=\"top\">Successful exit</td></tr>
  <tr><td valign=\"top\">info &gt; 0:</td>
      <td valign=\"top\" width=\"350\">If info = i then U[i,i] is exactly zero. The factorization
          has been completed, but the factor U is exactly
          singular, and division by zero will occur if it is used
          to solve a system of equations.</td></tr>
</table>
</blockquote>
<p>
The LU factorization is computed
with the LAPACK function &quot;dgetrf&quot;,
i.e., by Gaussian elemination using partial pivoting
with row interchanges. Vector &quot;pivots&quot; are the
pivot indices, i.e., for 1 &le; i &le; min(m,n), row i of
matrix A was interchanged with row pivots[i].
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real b1[3] = {10,22,12};
  Real b2[3] = { 7,13,10};
  Real    LU[3,3];
  Integer pivots[3];
  Real    x1[3];
  Real    x2[3];
<strong>algorithm</strong>
  (LU, pivots) := Matrices.LU(A);
  x1 := Matrices.LU_solve(LU, pivots, b1);  // x1 = {3,2,1}
  x2 := Matrices.LU_solve(LU, pivots, b2);  // x2 = {1,0,2}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.LU_solve\">Modelica.Math.Matrices.LU_solve</a>,
<a href=\"modelica://Modelica.Math.Matrices.solve\">Modelica.Math.Matrices.solve</a>
</p>
</html>"));
end LU;
