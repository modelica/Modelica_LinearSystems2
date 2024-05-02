within Modelica_LinearSystems2.Controllers.Internal;
function scaleFactor1 "Return scale factor for first order block"
  extends Modelica.Icons.Function;

  input Real n "(s+n)/(s+d)";
  input Real d "(s+n)/(s+d)";
  input Real small=100*Modelica.Constants.eps;
  output Real k "= d/n, if d,n are not zero, otherwise special cases";
algorithm
//  k :=(if abs(d) > small then abs(d) else 1)/(if abs(n) > small then abs(n) else 1);
  k := if abs(d) > small  and abs(n) > small then abs(d)/abs(n) else 1;
end scaleFactor1;
