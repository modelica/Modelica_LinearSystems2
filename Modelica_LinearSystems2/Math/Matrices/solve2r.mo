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
  output Real X[size(B, 1),size(B, 2)]=B "Matrix X such that X*op(A) = B";
  output Integer info;

protected
  Integer n=size(A, 1);
  Integer m=size(B, 2);
  Integer i;
  Integer j;
  Real LU[size(A, 1),size(A, 2)]=A;
  Integer k[size(A, 1)]= 1:n "Pivot vector";
  Real h[size(B, 1)];
//  Integer info;

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

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
X = Matrices.<strong>solve2r</strong>(A,B);
X = Matrices.<strong>solve2r</strong>(A, B, transA=false, isTriangular=false);
</pre></blockquote>

<h4>Description</h4>
<p>
This function call returns the
solution <strong>X</strong> of the linear system of equations
</p>
<blockquote>
<p>
<strong>X</strong>*op<strong>(A)</strong> = <strong>B</strong>
</p>
</blockquote>
<p>
with
</p>
<blockquote>
<p>
op<strong>(A)</strong> = transpose(<strong>(A)</strong>)  if   transA==true
op<strong>(A)</strong> = <strong>(A)</strong>  if   transA==false
</p>
</blockquote>
<p>
If matrix <strong>(A)</strong> is already lower triangular, the factorization is avoided if input \"isTriangular\" is set true.
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

  Real B[2,3]  = [10, 22, 12;
                  20, 44, 24];
  Real X[2,3];
<strong>algorithm</strong>
  X := Matrices.solve2r(A, B);  /* X = [-34.0, 17.2, 2.4;
                                        -68.0, 34.4, 4.8] */
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.LU_solve2\">Matrices.LU_solve2</a>
</p>
</html>"));
end solve2r;
