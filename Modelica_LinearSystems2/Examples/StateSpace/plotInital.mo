within Modelica_LinearSystems2.Examples.StateSpace;
function plotInital "Initial condition plot example"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);
protected
  Real x0[size(ss.A, 1)]=ones(size(ss.A, 1)) "Initial state vector";
algorithm
  Modelica_LinearSystems2.StateSpace.Plot.initialResponse(
     ss=ss, x0=x0);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the initial condition response of the system
StateSpace <em>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</em> to the initial condition <em>x0=[1;1]</em>.
</p>
</html>"));
end plotInital;
