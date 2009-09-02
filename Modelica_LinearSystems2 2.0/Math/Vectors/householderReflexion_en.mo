within Modelica_LinearSystems2.Math.Vectors;
function householderReflexion_en
  "reflect vector a to the unity vector e_n={0, ..., 0, 1}"
  import Modelica_LinearSystems2.Math.Vectors;

  input Real a[:];
  input Real u[size(a, 1)] "householder vector";
  output Real ra[size(u, 1)] "reflexion of a";

protected
  Integer n=size(u, 1);

algorithm
  ra := cat(
    1,
    fill(0, n - 1),
    {a[n] - 2*u[n]*(u*a)});

end householderReflexion_en;
