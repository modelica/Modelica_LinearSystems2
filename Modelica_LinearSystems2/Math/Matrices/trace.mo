within Modelica_LinearSystems2.Math.Matrices;
function trace "Sum of the diagonal elements of A"
  extends Modelica.Icons.Function;

  input Real A[:,size(A, 1)];
  output Real result;

protected
  Integer n=size(A, 1);
  Real r;

algorithm
  r := 0;
  if n > 0 then
    for i in 1:n loop
      r := r + A[i, i];
    end for;
  end if;
  result := r;
  annotation (
    obsolete = "Obsolete function - use Modelica.Math.Matrices.trace instead");
end trace;
