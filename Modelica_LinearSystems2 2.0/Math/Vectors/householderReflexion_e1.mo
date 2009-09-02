within Modelica_LinearSystems2.Math.Vectors;
function householderReflexion_e1
  "reflect vector a to the unity vector e_1={1,0,...,0}"
  import Modelica_LinearSystems2.Math.Vectors;

  input Real a[:];
  input Real u[size(a, 1)] "householder vector";
  output Real ra[size(u, 1)] "reflexion of a";

protected
  Integer n=size(u, 1);

algorithm
  ra := cat(
    1,
    {a[1] - 2*u[1]*(u*a)},
    fill(0, n - 1));

end householderReflexion_e1;
