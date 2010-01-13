within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_rank "Rank of a complex matrix (computed with singular values)"
  import Modelica_LinearSystems2;
  extends Modelica.Icons.Function;
  input Modelica_LinearSystems2.Math.Complex A[:,:] "Matrix";
  input Real eps=0
    "If eps > 0, the singular values are checked against eps; otherwise eps=max(size(A))*norm(A)*Modelica.Constants.eps is used";
  output Integer result "Rank of matrix A";

protected
  Integer n=min(size(A, 1), size(A, 2));
  Integer i=n;
  Real sigma[n];
  Real eps2;
algorithm
  result := 0;
  if n > 0 then
    sigma := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_singularValues(
                                                                    A);
    eps2 := if eps > 0 then eps else max(size(A))*sigma[1]*Modelica.Constants.eps;
    while i > 0 loop
      if sigma[i] > eps2 then
        result := i;
        i := 0;
      end if;
      i := i - 1;
    end while;
  end if;
  annotation (Documentation(info="<html>

</html>"));
end C_rank;
