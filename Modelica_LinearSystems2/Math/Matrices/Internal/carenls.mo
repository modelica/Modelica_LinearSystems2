within Modelica_LinearSystems2.Math.Matrices.Internal;
function carenls
  "Newton's method with exact line search for solving continuous algebraic riccati equation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Polynomial;

  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),:];
  input Real R[size(B, 2),size(B, 2)]=identity(size(B, 2));
  input Real Q[size(A, 1),size(A, 2)]=identity(size(A, 1));
  input Real X0[size(A, 1),size(A, 2)];
  input Real eps=Matrices.Internal.frobeniusNorm(A)*1e-9;

  output Real X[size(X0, 1),size(X0, 2)];
  output Real r;

protected
  Integer n=size(A, 1);
  Real G[size(A, 1),size(A, 2)]=B*Modelica.Math.Matrices.solve2(R, transpose(B));
  Real Xk[size(X, 1),size(X, 2)];
  Real Ak[size(A, 1),size(A, 2)];
  Real Rk[size(A, 1),size(A, 2)];
  Real Nk[size(A, 1),size(A, 2)];
  Real Vk[size(A, 1),size(A, 2)];
  Real tk;
  Integer k;
  Complex xc[2];
  Boolean stop;

algorithm
  if n > 1 then
    k := 0;
    stop := false;
    Xk := X0;
    while (not stop and k<10) loop
      k := k + 1;
      Ak := A - G*Xk;
      Rk := transpose(A)*Xk + Xk*A + Q - Xk*G*Xk;
      Nk := Matrices.lyapunov(Ak, -Rk);
      Vk := Nk*G*Nk;
      tk := Matrices.Internal.findLocal_tk(Rk, Vk);
      stop := eps > Matrices.Internal.frobeniusNorm(tk*Nk)/Matrices.Internal.frobeniusNorm(Xk);
      Xk := Xk + tk*Nk;
    end while;
    X := Xk;
    r := Matrices.Internal.frobeniusNorm(X*A + transpose(A)*X - X*G*X + Q);

  elseif n == 1 then
    xc := Polynomial.roots(Polynomial({-G[1, 1],2*A[1, 1],Q[1, 1]}));
    X := matrix(-abs(xc[1].re));
    r := 0;
  else
    X := fill(0, 0, 0);
    r := 0;
  end if;

  annotation (Documentation(revisions="<html>
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
end carenls;
