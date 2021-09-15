within Modelica_LinearSystems2.Math.Matrices;
function solve2
  "Solve real system of linear equations A*X=B with a B matrix (Gaussian elemination with partial pivoting)"

  extends Modelica.Icons.Function;
  input Real A[:,size(A, 1)] "Matrix A of A*X = B";
  input Real B[size(A, 1),:] "Matrix B of A*X = B";
  output Real X[size(B, 1),size(B, 2)] "Matrix X such that A*X = B";

protected
  Integer info;
algorithm
  if size(A, 1) > 0 then
    (X,info) := Modelica.Math.Matrices.LAPACK.dgesv(A, B);
    assert(info == 0, "Solving a linear system of equations with function
\"Matrices.solve2\" is not possible, because the system has either
no or infinitely many solutions (A is singular).");
  else
    X := fill(
      0,
      size(B, 1),
      size(B, 2));
  end if;
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.solve2 instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<strong>solve2</strong>(A,b);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
solution <strong>X</strong> of the linear system of equations
</p>
<blockquote>
<p>
<strong>A</strong>*<strong>X</strong> = <strong>B</strong>
</p>
</blockquote>
<p>
If a unique solution <strong>X</strong> does not exist (since <strong>A</strong> is singular),
an exception is raised.
</p>

<h4>Note</h4>
<p>
The solution is computed with the LAPACK function \"dgesv\",
i.e., by Gaussian elemination with partial pivoting.
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real B[3,2] = [10, 20;
                 22, 44;
                 12, 24];
  Real X[3,2];
<strong>algorithm</strong>
  (LU, pivots) := Matrices.LU(A);
  X := Matrices.solve2(A, B1);  /* X = [3, 6;
                                        2, 4;
                                        1, 2] */
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.LU_solve2\">Matrices.LU_solve2</a>
</p>
</html>"));
end solve2;
