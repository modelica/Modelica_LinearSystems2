within Modelica_LinearSystems2.Math.Matrices;
function sylvester
  "Solution of continuous-time Sylvester equation A*X + X*B = C"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:] "Matrix A";
  input Real B[:,:] "Matrix B";
  input Real C[size(A, 1),size(B, 2)] "Matrix C";
  input Boolean aIsSchur=false "True if A has already real Schur form";
  input Boolean bIsSchur=false "True if B has already real Schur form";
  output Real X[size(A, 1),size(B, 2)] "Solution of Sylvester equation";

protected
  Integer n=size(A, 1);
  Integer m=size(B, 1);
  Real S[size(A, 1),size(A, 2)];
  Real T[size(B, 1),size(B, 2)];
  Real U[size(A, 1),size(A, 1)];
  Real V[size(B, 1),size(B, 1)];
  Real Chat[size(C, 1),size(C, 2)];
  Real scale;
  Integer info;

algorithm
  if n > 1 and m > 1 then
    if aIsSchur then
      S := A;
      U := identity(n);
    else
      (S,U) := Matrices.rsf2(A);
    end if;
    if bIsSchur then
      T := B;
      V := identity(m);
    else
      (T,V) := Matrices.rsf2(B);
    end if;

    Chat := if aIsSchur and bIsSchur then C else if aIsSchur then C*V else if
      bIsSchur then transpose(U)*C else transpose(U)*C*V;
    (X,scale,info) := Matrices.LAPACK.dtrsyl(S, T, Chat);
    assert(info == 0, "Solving of Sylvester equation with Matrices.sylvester was not successfull.\n
                    The value of info is " + String(info) + ", but should be zero. A value unequal to zero means:\n
            < 0: if INFO = -i, the i-th argument had an illegal value\n
            = 1: A and B have common or very close eigenvalues; perturbed
                 values were used to solve the equation (but the matrices
                 A and B are unchanged).");
    X := if aIsSchur and bIsSchur then scale*X else if aIsSchur then scale*X*
      transpose(V) else if bIsSchur then scale*U*X else scale*U*X*transpose(V);
  else
    X := fill(0, n, m);
  end if;

  annotation (Documentation(info="<html>
<p>
This function computes the solution <b>X</b> of the continuous-time Sylvester equation
</p>
<blockquote><pre>
<b>A</b>*<b>X</b> + <b>X</b>*<b>B</b> = <b>C</b>
</pre></blockquote>
<p>
using the Schur method for Sylvester equations proposed by Bartels and Stewart [1].
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Bartels, R.H. and Stewart G.W. (1972):</dt>
<dd> <b>Algorithm 432: Solution of the matrix equation AX + XB = C</b>.
     Comm. ACM., Vol. 15, pp. 820-826.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end sylvester;
