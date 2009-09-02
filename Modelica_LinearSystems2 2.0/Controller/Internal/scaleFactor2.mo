within Modelica_LinearSystems2.Controller.Internal;
function scaleFactor2 "Return scale factor for second order block"

  input Real n1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
  input Real n2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
  input Real d1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
  input Real d2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
  input Real small=100*Modelica.Constants.eps;
  output Real k "= d2/n2, if d2,n2 are not zero, otherwise special cases";
algorithm
  k :=(if abs(d2) > small then abs(d2) else (if abs(d1) > small then abs(d1) else 1)) /
      (if abs(n2) > small then abs(n2) else (if abs(n1) > small then abs(n1) else 1));
end scaleFactor2;
