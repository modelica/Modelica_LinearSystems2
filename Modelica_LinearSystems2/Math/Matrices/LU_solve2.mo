within Modelica_LinearSystems2.Math.Matrices;
function LU_solve2
  "Solve real system of linear equations P*L*U*X=B with a B vector and an LU decomposition (from LU(..))"

  input Real LU[:,size(LU, 1)]
    "L,U factors of Matrices.LU(..) for a square matrix";
  input Integer pivots[size(LU, 1)] "Pivots indices of Matrices.LU(..)";
  input Real B[size(LU, 1),:] "Right hand side matrix of P*L*U*X=B";
  output Real X[size(B, 1),size(B, 2)] "Solution matrix such that P*L*U*X = B";

algorithm
  for i in 1:size(LU, 1) loop
    assert(LU[i, i] <> 0, "Solving a linear system of equations with function
\"Matrices.LU_solve\" is not possible, since the LU decomposition
is singular, i.e., no unique solution exists.");
  end for;
  if size(LU, 1) > 0 then
    X := LAPACK.dgetrs(
      LU,
      pivots,
      B);
  else
    X := fill(
      0,
      size(B, 1),
      size(B, 2));
  end if;
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.LU_solve2 instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<b>LU_solve2</b>(LU, pivots, B);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
solution <b>x</b> of the linear systems of equations
</p>
<blockquote>
<p>
<b>P</b>*<b>L</b>*<b>U</b>*<b>X</b> = <b>B</b>;
</p>
</blockquote>
<p>
where <b>P</b> is a permutation matrix (implicitely
defined by vector <code>pivots</code>),
<b>L</b> is a lower triangular matrix with unit
diagonal elements (lower trapezoidal if m &gt; n), and
<b>U</b> is an upper triangular matrix (upper trapezoidal if m &lt; n).
The matrices of this decomposition are computed with function
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a> that
returns arguments <code>LU</code> and <code>pivots</code>
used as input arguments of <code>Matrices.LU_solve</code>.
With <code>Matrices.LU</code> and <code>Matrices.LU_solve</code>
it is possible to efficiently solve linear systems
with different right hand side <b>matrices</b>. If a linear system of equations with
just one right hand side matrix shall be solved, it is
more convenient to just use the function
<a href=\"modelica://Modelica.Math.Matrices.solve2\">Matrices.solve2</a>.
</p>
<p>
If a unique solution <b>X</b> does not exist (since the
LU decomposition is singular), an exception is raised.
</p>
<p>
The LU factorization is computed
with the LAPACK function \"dgetrf\",
i.e., by Gaussian elemination using partial pivoting
with row interchanges. Vector \"pivots\" are the
pivot indices, i.e., for 1 &le; i &le; min(m,n), row i of
matrix A was interchanged with row pivots[i].
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real B1[3] = [10, 20;
                22, 44;
                12, 24];
  Real B2[3] = [ 7, 14;
                13, 26;
                10, 20];
  Real    LU[3,3];
  Integer pivots[3];
  Real    X1[3,2];
  Real    X2[3,2];
<b>algorithm</b>
  (LU, pivots) := Matrices.LU(A);
  X1 := Matrices.LU_solve2(LU, pivots, B1);  /* X1 = [3, 6;
                                                      2, 4;
                                                      1, 2] */
  X2 := Matrices.LU_solve2(LU, pivots, B2);  /* X2 = [1, 2;
                                                      0, 0;
                                                      2, 4] */
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.solve2\">Matrices.solve2</a>
</p>
</html>"));
end LU_solve2;
