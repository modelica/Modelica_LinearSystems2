within Modelica_LinearSystems2.Math.Matrices;
function dlyapunov
  "Solution of continuous-time Lyapunov equation A'X*A - X = C"
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Math.Matrices.solve;

  input Real A[:,size(A, 1)];
  input Real C[size(A, 1),size(A, 2)];
  input Real eps=Modelica.Math.Matrices.norm(A,1)*10*Modelica.Constants.eps;

protected
  Integer n=size(A, 1);
  Real R[size(A, 1),size(A, 2)] "rsf of A', i.e. R=U'A'U";
  Real U[size(A, 1),size(A, 2)] "transformation matrix U for R=U'A'U";
  Real C2[size(A, 1),size(A, 2)];
  Real R22[size(A, 1),size(A, 2)];
  Real R11[size(A, 1),size(A, 2)];
  Integer k;

//  Real Y[size(A, 1),size(A, 2)];
  Real bk1[size(A, 1)];
  Real bk[size(A, 1)];
  Real yk[size(A, 1)];
  Boolean crit;

public
  output Real X[size(A, 1),size(A, 2)] "solution of the Lyapunov equation";

algorithm
     X := zeros(n, n);
//    Y := zeros(n,n);
    k := n;
  if n > 1 then
    (R,U) := Matrices.rsf(transpose(A));
    C2 := transpose(U)*C*U;

  while k > 0 loop

//      bk := C2[:,k] - R*Y[:,k+1:n]*R[k,k+1:n];
      bk := C2[:,k] - R*X[:,k+1:n]*R[k,k+1:n];
      crit :=  if k>1 then abs(R[k,k-1])<eps else false;
   if (k==1 or crit) then
      R22 := R[k,k]*R;
      for i in 1:n loop
        R22[i,i] := R22[i,i]-1.0;
      end for;
//      Y[:,k] := solve(R22,bk);
      X[:,k] := solve(R22,bk);
      k:=k-1;
   else
//       bk1:=C2[:,k-1]-R*Y[:,k+1:n]*R[k-1,k+1:n];
       bk1:=C2[:,k-1]-R*X[:,k+1:n]*R[k-1,k+1:n];
       R11 := R[k-1,k-1]*R;
       R22 :=  R[k,k]*R;
       for i in 1:n loop
         R11[i,i] := R11[i,i] -1.0;
         R22[i,i] := R22[i,i] -1.0;
       end for;
       yk :=solve([R11,  R[k-1,k]*R; R[k,k-1]*R,  R22], cat(1,bk1,bk));
//       Y[:,k-1]:=yk[1:n];
//       Y[:,k]:=yk[n + 1:2*n];
       X[:,k-1]:=yk[1:n];
       X[:,k]:=yk[n + 1:2*n];
       k:=k-2;
   end if;
   end while;

// transform X corresponding to the original form
//    X := U*Y*transpose(U);
    X := U*X*transpose(U);

  elseif n == 1 then
    X[1, 1] := C[1, 1]/(A[1, 1]*A[1, 1]-1);
  else
    X := fill(0, 0, 0);
  end if;

  annotation (Documentation(info="<html>
<p>
Function <b>laypunov</b> computes the solution <b>X</b> of the continuous-time
Lyapunov equation
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
</html>"));
end dlyapunov;
