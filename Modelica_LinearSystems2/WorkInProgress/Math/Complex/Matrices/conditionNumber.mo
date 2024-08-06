within Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices;
encapsulated function conditionNumber
  "Calculate the condition number norm(A)*norm(inv(A))"
  extends Modelica.Icons.Function;
  import Complex;
  import Modelica_LinearSystems2;
  import Modelica;

  input Complex A[:,:] "Input matrix";
  input Real p(min=1) = 2
    "Type of p-norm (only allowed: 1, 2 or Modelica.Constants.inf)";
  output Real result=0.0 "p-norm of matrix A";
protected
  Real eps=1e-20;
  Real s[size(A, 1)] "singular values";

algorithm
  if p == 2 then
    s := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_singularValues(A);
    if min(s) < eps then
      result := -1e100;
    else
      result := max(s)/min(s);
    end if;
  else
    result := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.norm(
                                    A, p)*Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.norm(
                                                                Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.inv(
                                                                                     A), p);
  end if;
end conditionNumber;
