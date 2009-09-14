within Modelica_LinearSystems2.Examples.StateSpace;
function designAssignPolesSISO
  "Example for pole placing using Ackermann's method"
  extends Modelica.Icons.Function;
  output Real k[2] "Gain vector";
  annotation (interactive=true, Documentation(info="<html>
<p>
Computes the gain vector k for the state space system
<pre> 
sc = StateSpace(A=[-1,1;0,-2],B=[0, 1],C=[1,0; 0, 1],D=[0; 0])
</pre>
such that for the state feedback 
<pre>u = -k*y = -k*x</pre> the closed-loop
poles are placed at 
<pre>p = {-3,-4}.</pre>
</html>"));

  input Modelica_LinearSystems2.StateSpace sc=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[0; 1],
      C=[1,0; 0,1],
      D=[0; 0]);
protected
      Modelica_LinearSystems2.Math.Complex p[2]={
      Modelica_LinearSystems2.Math.Complex(-3,0),
      Modelica_LinearSystems2.Math.Complex(-4,0)};
algorithm
  k := Modelica_LinearSystems2.StateSpace.Design.assignPolesSI(sc, p);
end designAssignPolesSISO;
