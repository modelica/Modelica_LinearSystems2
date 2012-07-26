within Modelica_LinearSystems2.Math.Matrices;
function lyapunov
  "Solution of continuous-time Lyapunov equation X*A + A'*X = C"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,size(A, 1)];
  input Real C[size(A, 1),size(A, 2)];
  input Real eps=Modelica.Math.Matrices.norm(A,1)*10*Modelica.Constants.eps;

protected
  Integer n=size(A, 1);
  Real R[size(A, 1),size(A, 2)] "rsf of A', i.e. R=U'A'U";
  Real U[size(A, 1),size(A, 2)] "transformation matrix U for R=U'A'U";
  Real C2[size(A, 1),size(A, 2)];
  Real R11[size(A, 1),size(A, 2)];
  Real R22[size(A, 1),size(A, 2)];
  Real R12[size(A, 1),size(A, 2)];
  Real R21[size(A, 1),size(A, 2)];
  Real R2[2*size(A, 1),2*size(A, 2)];
  Real I[size(A, 1),size(A, 1)]=identity(size(A, 1));
  Real x[2*size(A,1)];
  Real c[2*size(A,1)];
  Real CC[size(A,1),2];
  Integer k;

public
  output Real X[size(A, 1),size(A, 2)] "solution of the Lyapunov equation";

algorithm
  if n > 1 then
    (R,U) := Matrices.rsf(transpose(A));
    C2 := transpose(U)*C*U;
    X := zeros(n, n);

// Calculate the last 1 or 2 columns of X
    R22 := R;
    for i1 in 1:n loop
      R22[i1, i1] := R[i1, i1] + R[n, n];
    end for;
    if abs(R[n, n - 1]) < eps then
      X[:, n] := Matrices.solve(R22, C2[:, n]);
      k := n - 1;
    else
      R11 := R;
      R12 := zeros(n, n);
      R21 := zeros(n, n);
      for i1 in 1:n loop
        R11[i1, i1] := R[i1, i1] + R[n - 1, n - 1];
        R12[i1, i1] := R[n - 1, n];
        R21[i1, i1] := R[n, n - 1];
      end for;

// solve 2nx2n equation for 2x2 Schur bump with Kronecker product and vec operator approach
      R2 := [R11,R12; R21,R22];
      c := cat(1, C2[:, n - 1], C2[:, n]);
      x := Matrices.solve(R2, c);
      X[:, n - 1] := x[1:n];
      X[:, n] := x[n + 1:2*n];
      k := n - 2;
    end if;

// Calculate the rest of X

    while k > 1 loop
      R22 := R;
      for i1 in 1:n loop
        R22[i1, i1] := R[i1, i1] + R[k, k];
      end for;
      if abs(R[k, k - 1]) < eps then
        //real eigenvalue
//        X[:, k] := Matrices.solve(R22, C2[:, k] - vector(X[:, k + 1:n]*transpose(matrix(R[k, k + 1:n]))));
        X[:, k] := Matrices.solve(R22, C2[:, k] - vector(X[:, k + 1:n]*matrix(R[k, k + 1:n])));
        k := k - 1;
      else
       // conjugated complex eigenvalues
        R11 := R;
        R12 := zeros(n, n);
        R21 := zeros(n, n);
        for i1 in 1:n loop
          R11[i1, i1] := R[i1, i1] + R[k - 1, k - 1];
          R12[i1, i1] := R[k - 1, k];
          R21[i1, i1] := R[k, k - 1];
        end for;
        R2 := [R11,R12; R21,R22];
        CC := C2[:,k-1:k] - X[:,k+1:n]*transpose(R[k-1:k,k+1:n]);
        c := cat(1, CC[:, 1], CC[:, 2]);
        x := Matrices.solve(R2, c);
        X[:, k - 1] := x[1:n];
        X[:, k] := x[n + 1:2*n];

        k := k - 2;
      end if;
    end while;// k=1 or k=0

// if k=1 the first column (if there exist a real eigenvalue) has to be calculated separately
    if k == 1 then
      R22 := R;
      for i1 in 1:n loop
        R22[i1, i1] := R[i1, i1] + R[1, 1];
      end for;
      X[:, 1] := Matrices.solve(R22, C2[:, 1] - vector(X[:, 2:n]*matrix(R[1, 2:n])));
    end if;

// transform X corresponding to the original form
    X := U*X*transpose(U);

  elseif n == 1 then
    X[1, 1] := C[1, 1]/(2*A[1, 1]);
  else
    X := fill(0, 0, 0);
  end if;

  annotation (Documentation(info="<html>
<p>
Function <b>laypunov</b> computes the solution <b>X</b> of the continuous-time Lyapunov equation
</p>
<blockquote><pre>
<b>X</b><b>A</b> + <b>A</b>'*<b>X</b> = <b>C</b>.
</pre></blockquote>
<p>
using the Schur method for Lyapunov equations proposed by Bartels and Stewart [1].
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
</html>
"));
end lyapunov;
