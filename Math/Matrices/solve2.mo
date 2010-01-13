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
  annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
Matrices.<b>solve2</b>(A,b);
</pre></blockquote>
<h4>Description</h4>
<p>
This function call returns the
solution <b>X</b> of the linear system of equations
</p>
<blockquote>
<p>
<b>A</b>*<b>X</b> = <b>B</b>
</p>
</blockquote>
<p>
If a unique solution <b>X</b> does not exist (since <b>A</b> is singular),
an exception is raised.
</p>
<p>
Note, the solution is computed with the LAPACK function \"dgesv\",
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
<b>algorithm</b>
  (LU, pivots) := Matrices.LU(A);
  X := Matrices.solve2(A, B1);  /* X = [3, 6;
                                        2, 4;
                                        1, 2] */
</pre></blockquote>
 
<h4>See also</h4>
<a href=\"Modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"Modelica://Modelica.Math.Matrices.LU_solve2\">Matrices.LU_solve2</a>
</HTML>"));
end solve2;
