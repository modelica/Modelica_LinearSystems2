within Modelica_LinearSystems2.Math.Matrices.Internal;
function findLocal_tk
  "Find a local minimizer tk to define the length of the step tk*Nk in carenls"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2.Math.Complex;

  input Real Rk[:,size(Rk, 2)];
  input Real G[size(Rk, 1),size(Rk, 2)];
  input Real Nk[size(Rk, 1),size(Rk, 2)];

  output Real tk;

protected
  Real Vk[size(Rk, 1),size(Rk, 2)];
  Real alpha_k;
  Real beta_k;
  Real gamma_k;
  Complex p[3];
  Boolean h;

algorithm
  Vk := Nk*G*Nk;
  alpha_k := Matrices.trace(Rk*Rk);
  beta_k := Matrices.trace(Rk*Vk);
  gamma_k := Matrices.trace(Vk*Vk);

  if gamma_k > Modelica.Constants.eps then
    p := Polynomial.roots(Polynomial({4*gamma_k,6*beta_k,2*(alpha_k - 2*beta_k),
      -2*alpha_k}));
    h := false;
    for i1 in 1:3 loop
      if (abs(Complex.imag(p[i1])) < Modelica.Constants.eps) then
        if (abs(Complex.real(p[i1]) - 1) <= 1) then
          tk := Complex.real(p[i1]);
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

end findLocal_tk;
