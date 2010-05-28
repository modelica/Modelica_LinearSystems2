within Modelica_LinearSystems2.Math.Matrices;
function solve2r
  "Solve real system of linear equations X*op(A)=B with a B matrix (Gaussian elemination with partial pivoting)"

  extends Modelica.Icons.Function;
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real A[:,size(A,1)] "Matrix A of X*op(A) = B";
  input Real B[:,size(A,1)] "Matrix B of X*op(A) = B";
  input Boolean transA=false "True if op(A)=A', false if op(A)=A";
  input Boolean isTriangular=false "True if the A is already lower triangular";
  output Real X[size(B, 1),size(A, 1)]=B "Matrix X such that X*op(A) = B";

protected
  Integer n=size(A, 1);
  Integer m=size(B, 1);
  Integer i;
  Integer j;
  Real LU[size(A, 1),size(A, 2)]=A;
  Integer k[size(A, 1)] "Pivot vector";
  Real h[size(B, 1)];
  Integer info;

algorithm
  if not isTriangular then
    (LU,k,info) := Modelica.Math.Matrices.LAPACK.dgetrf(A);
    assert(info == 0, "LU factoriztion failed in \"solve2r\"");
  else
    LU := A;
  end if;

  if transA then
    for i in 1:n loop
      j := k[i];
      if (j <> i) then
        h := X[:, j];
        X[:, j] := X[:, i];
        X[:, i] := h;
      end if;
    end for;
    X := LAPACK.dtrsm(LU, X, 1, true, false, true, true);
    X := LAPACK.dtrsm(LU, X, 1, true, true, true, false);
  else
    X := LAPACK.dtrsm(LU, X, 1, true, true, false, false);
    X := LAPACK.dtrsm(LU, X, 1, true, false, false,
      true);
    for i in n:-1:1 loop
      j := k[i];
      if (j <> i) then
        h := X[:, j];
        X[:, j] := X[:, i];
        X[:, i] := h;
      end if;
    end for;
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
end solve2r;
