within Modelica_LinearSystems2.Math.Vectors;
function householderReflexion
  "Reflect vector a on a plane with orthogonal vector u"
  extends Modelica.Icons.Function;

  input Real a[:];
  input Real u[size(a, 1)] "Householder vector";
  output Real ra[size(u, 1)] "Reflexion of a";

protected
  Real norm_a=Modelica.Math.Vectors.length(a);
  Real h=2*u*a;

algorithm
  ra := a - h*u;

// this function is mainly used in the fromStateSpace transformations.
// In this context the calculation of invariant zeros is very susceptible with subject to
// elements of the output vector are zero or not. Therefore, values close to zero are set to zero.
  for i in 1:size(ra, 1) loop
    ra[i] := if abs(ra[i]) >= norm_a*1e-12 then ra[i] else 0;
  end for;

  annotation (Documentation(info="<html>
</html>"));
end householderReflexion;
