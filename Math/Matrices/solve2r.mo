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

  annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
  X = Matrices.<b>solve2r</b>(A,B);
  X = Matrices.<b>solve2r</b>(A, B, transA=false, isTriangular=false);
</pre></blockquote>
<h4>Description</h4>
<p>
This function call returns the
solution <b>X</b> of the linear system of equations
</p>
<blockquote>
<p>
<b>X</b>*op<b>(A)</b> = <b>B</b>
</p>
</blockquote>
<p>
with 
</p>
<blockquote>
<p>
  op<b>(A)</b> = transpose(<b>(A)</b>)  if   transA==true
  op<b>(A)</b> = <b>(A)</b>  if   transA==false
</p>
</blockquote>
<p>
If matrix <b>(A)</b> is already lower triangular, the factorization is avoided if input \"isTriangular\" is set true.
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
  
  Real B[2,3]  = [10, 22, 12;
                  20, 44, 24];
  Real X[2,3];
<b>algorithm</b>
  X := Matrices.solve2r(A, B);  /* X = [-34.0, 17.2, 2.4;
                                        -68.0, 34.4, 4.8] */
</pre></blockquote>

<h4>See also</h4>
<a href=\"modelica://Modelica.Math.Matrices.LU\">Matrices.LU</a>,
<a href=\"modelica://Modelica.Math.Matrices.LU_solve2\">Matrices.LU_solve2</a>
</HTML>"));
end solve2r;
