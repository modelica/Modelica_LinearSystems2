within Modelica_LinearSystems2.Math.Matrices;
function LU_solve
  "Solve real system of linear equations P*L*U*x=b with a b vector and an LU decomposition (from LU(..))"

  extends Modelica.Icons.Function;
  input Real LU[:,size(LU, 1)]
    "L,U factors of Matrices.LU(..) for a square matrix";
  input Integer pivots[size(LU, 1)] "Pivots indices of Matrices.LU(..)";
  input Real b[size(LU, 1)] "Right hand side vector of P*L*U*x=b";
  output Real x[size(b, 1)] "Solution vector such that P*L*U*x = b";

algorithm
  for i in 1:size(LU, 1) loop
    assert(LU[i, i] <> 0, "Solving a linear system of equations with function
\"Matrices.LU_solve\" is not possible, since the LU decomposition
is singular, i.e., no unique solution exists.");
  end for;
  if size(LU, 1) > 0 then
    x := Modelica.Math.Matrices.LAPACK.dgetrs_vec(
      LU,
      pivots,
      b);
  else
    x := fill(0, 0);
  end if;
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.LU_solve instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<strong>LU_solve</strong>(LU, pivots, b);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
solution <strong>x</strong> of the linear systems of equations
</p>
<blockquote>
<p>
<strong>P</strong>*<strong>L</strong>*<strong>U</strong>*<strong>x</strong> = <strong>b</strong>;
</p>
</blockquote>
<p>
where <strong>P</strong> is a permutation matrix (implicitely
defined by vector <code>pivots</code>),
<strong>L</strong> is a lower triangular matrix with unit
diagonal elements (lower trapezoidal if m &gt; n), and
<strong>U</strong> is an upper triangular matrix (upper trapezoidal if m &lt; n).
The matrices of this decomposition are computed with function
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a> that
returns arguments <code>LU</code> and <code>pivots</code>
used as input arguments of <code>Matrices.LU_solve</code>.
With <code>Matrices.LU</code> and <code>Matrices.LU_solve</code>
it is possible to efficiently solve linear systems
with different right hand side vectors. If a linear system of equations with
just one right hand side vector shall be solved, it is
more convenient to just use the function
<a href=\"modelica://Modelica.Math.Matrices.solve\">Matrices.solve</a>.
</p>
<p>
If a unique solution <strong>x</strong> does not exist (since the
LU decomposition is singular), an exception is raised.
</p>
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
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.solve\">Matrices.solve</a>
</p>
</html>"));
end LU_solve;
