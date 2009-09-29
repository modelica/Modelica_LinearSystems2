within Modelica_LinearSystems2.Math.Matrices.Internal;
function conditionNumber "Calculate the condition number norm(A)*norm(inv(A))"
  extends Modelica.Icons.Function;

  input Real A[:,:] "Input matrix";
  input Real p(min=1) = 2
    "Type of p-norm (only allowed: 1, 2 or Modelica.Constants.inf)";
  output Real result=0.0 "p-norm of matrix A";

protected
  Real eps=1e-20;
  Real s[size(A, 1)] "singular values";

algorithm
  if p == 2 then
    s := Modelica.Math.Matrices.singularValues(A);
    if min(s) < eps then
result := -1e100;
    else
result := max(s)/min(s);
    end if;
  else
    result := Modelica.Math.Matrices.norm(A, p)*Modelica.Math.Matrices.norm(
      Modelica.Math.Matrices.inv(A), p);
  end if;
end conditionNumber;
