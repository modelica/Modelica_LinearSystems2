within Modelica_LinearSystems2.Math.Matrices;
function solve
  "Solve real system of linear equations A*x=b with a b vector (Gaussian elemination with partial pivoting)"

  extends Modelica.Icons.Function;
  input Real A[:,size(A, 1)] "Matrix A of A*x = b";
  input Real b[size(A, 1)] "Vector b of A*x = b";
  output Real x[size(b, 1)] "Vector x such that A*x = b";

protected
  Integer info;
algorithm
  if size(A, 1) > 0 then
    (x,info) := Modelica.Math.Matrices.LAPACK.dgesv_vec(A, b);
    assert(info == 0, "Solving a linear system of equations with function
\"Matrices.solve\" is not possible, because the system has either
no or infinitely many solutions (A is singular).");
  else
    x := fill(0, 0);
  end if;
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.solve instead",
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<strong>solve</strong>(A,b);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
solution <strong>x</strong> of the linear system of equations
</p>
<blockquote>
<strong>A</strong>*<strong>x</strong> = <strong>b</strong>
</blockquote>
<p>
If a unique solution <strong>x</strong> does not exist (since <strong>A</strong> is singular),
an exception is raised.
</p>

<h4>Note</h4>
<p>
The solution is computed with the LAPACK function &quot;dgesv&quot;,
i.e., by Gaussian elemination with partial pivoting.
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1,2,3;
                 3,4,5;
                 2,1,4];
  Real b[3] = {10,22,12};
  Real x[3];
<strong>algorithm</strong>
  x := Matrices.solve(A,b);  // x = {3,2,1}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.LU\">Modelica.Math.Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.LU_solve\">Modelica.Math.Matrices.LU_solve</a>
</p>
</html>"));
end solve;
