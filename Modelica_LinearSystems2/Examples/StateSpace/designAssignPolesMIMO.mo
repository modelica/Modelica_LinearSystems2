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
  input Complex p[:]={Complex(-3, 0),Complex(-4, 0)};

protected
  Complex newPoles[:];

algorithm
  (K, S, newPoles) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, p);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("newPoles", newPoles);
  newPoles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A - ss.B*K);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("newPoles", newPoles);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the gain vector k for the state space system
</p>
<blockquote><pre>
ss = StateSpace(A=[-1,1;0,-2],B=[0, 1],C=[1,0; 0, 1],D=[0; 0])
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
end designAssignPolesMIMO;
