within Modelica_LinearSystems2.Examples.StateSpace;
function designAssignPolesMIMO "Example for pole placing"
  extends Modelica.Icons.Function;
  output Real K[:,:] "Gain vector";
  output Real S[:,:];

  input Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[0; 1],
      C=[1,0; 0,1],
      D=[0; 0]);
      input Modelica_LinearSystems2.Math.Complex p[:]={
      Modelica_LinearSystems2.Math.Complex(-3,0),
      Modelica_LinearSystems2.Math.Complex(-4,0)};

protected
      Modelica_LinearSystems2.Math.Complex newPoles[:];

algorithm
  (K, S, newPoles) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, p);
  Modelica_LinearSystems2.Math.Complex.Vectors.print("newPoles", newPoles);
  newPoles := Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A-ss.B*K);
  Modelica_LinearSystems2.Math.Complex.Vectors.print("newPoles", newPoles);

  annotation (interactive=true, Documentation(info="<html>
<p>
Computes the gain vector k for the state space system
<pre> 
ss = StateSpace(A=[-1,1;0,-2],B=[0, 1],C=[1,0; 0, 1],D=[0; 0])
</pre>
such that for the state feedback 
<pre>u = -k*y = -k*x</pre> the closed-loop
poles are placed at 
<pre>p = {-3,-4}.</pre>
</html>"));
end designAssignPolesMIMO;
