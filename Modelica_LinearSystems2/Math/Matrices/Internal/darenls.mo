within Modelica_LinearSystems2.Math.Matrices.Internal;
function darenls
  "Newton's method with exact line search for solving continuous algebraic riccati equation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;

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
  Real Xk[size(X, 1),size(X, 2)];
  Real Ak[size(A, 1),size(A, 2)];
  Real Rk[size(A, 1),size(A, 2)];
  Real Nk[size(A, 1),size(A, 2)];
  Real Hk[size(B, 2),size(B, 1)];
  Real Sk[size(B, 1),size(B, 1)];
  Real Vk[size(A, 1),size(A, 2)];
  Real tk;
  Integer k;
  Real AT[size(A, 2),size(A, 2)]=transpose(A);
  Real BT[size(B, 2),size(B, 1)]=transpose(B);

  Boolean stop;

algorithm
  if n > 0 then
    k := 0;
    stop := false;
    Xk := X0;
    while (not stop and k<10) loop
      k := k + 1;
      Hk := Modelica.Math.Matrices.solve2(R+BT*Xk*B,BT);
      Ak := A-B*Hk*Xk*A;

      Rk:=AT*Xk*A - Xk + Q - AT*Xk*B*Hk*Xk*A;

      Nk := Matrices.dlyapunov(Ak, -Rk);
//  Matrices.printMatrix(Nk,6,"Nk");

      Sk := B*Hk;
      Vk :=transpose(Ak)*Nk*Sk*Nk*Ak;
//   Matrices.printMatrix(Vk,6,"Vk");
      tk := Matrices.Internal.findLocal_tk(Rk, Vk);
      stop := eps > Matrices.Internal.frobeniusNorm(tk*Nk)/Matrices.Internal.frobeniusNorm(Xk);
      Xk := Xk + tk*Nk;

    end while;
    X := Xk;
    r := Matrices.Internal.frobeniusNorm(AT*X*A - X +Q - AT*X*B*Modelica.Math.Matrices.solve2(R+BT*X*B,BT)*X*A);
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
end darenls;
