within Modelica_LinearSystems2.Math.Matrices.Internal;
function findLocal_tk
  "Find a local minimizer tk to define the length of the step tk*Nk in carenls or darenls"
  extends Modelica.Icons.Function;

  import Modelica.ComplexMath;
  import Modelica_LinearSystems2.Math.Polynomial;

  input Real Rk[:,size(Rk, 2)];
  input Real Vk[size(Rk, 1),size(Rk, 2)];
//   input Real G[size(Rk, 1),size(Rk, 2)];
//   input Real Nk[size(Rk, 1),size(Rk, 2)];

  output Real tk;

//  Real Vk[size(Rk, 1),size(Rk, 2)];
protected
  Real alpha_k;
  Real beta_k;
  Real gamma_k;
  Complex p[3];
  Boolean h;

algorithm
//  Vk := Nk*G*Nk;
  alpha_k := Modelica.Math.Matrices.trace(Rk*Rk);
  beta_k := Modelica.Math.Matrices.trace(Rk*Vk);
  gamma_k := Modelica.Math.Matrices.trace(Vk*Vk);

  if gamma_k > Modelica.Constants.eps then
    p := Polynomial.roots(Polynomial({4*gamma_k,6*beta_k,2*(alpha_k - 2*beta_k),
      -2*alpha_k}));
    h := false;
    for i1 in 1:3 loop
      if (abs(ComplexMath.imag(p[i1])) < Modelica.Constants.eps) then
        if (abs(ComplexMath.real(p[i1]) - 1) <= 1) then
          tk := ComplexMath.real(p[i1]);
          h := true;
        end if;
      end if;
    end for;
    if not h then
      tk := 1;
    end if;

  else
    tk := 1;
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
end findLocal_tk;
