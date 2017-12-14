within Modelica_LinearSystems2.Examples.StateSpace;
function designAssignPolesSISO "Example for pole placing using Ackermann's method"
  extends Modelica.Icons.Function;

  output Real k[2] "Gain vector";

  input Modelica_LinearSystems2.StateSpace sc=Modelica_LinearSystems2.StateSpace(
    A=[-1,1; 0,-2],
    B=[0; 1],
    C=[1,0; 0,1],
    D=[0; 0]);
protected
  Complex p[2]={Complex(-3, 0),Complex(-4, 0)};
algorithm
  k := Modelica_LinearSystems2.StateSpace.Design.assignPolesSI(sc, p);

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Computes the gain vector k for the state space system
</p>
<blockquote><pre>
sc = StateSpace(A=[-1,1;0,-2],B=[0, 1],C=[1,0; 0, 1],D=[0; 0])
</pre></blockquote>
<p>
such that for the state feedback
</p>
<blockquote><pre>
u = -k*y = -k*x
</pre></blockquote>
<p>
the closed-loop poles are placed at
</p>
<blockquote><pre>
p = {-3,-4}.
</pre></blockquote>
</html>"));
end designAssignPolesSISO;
