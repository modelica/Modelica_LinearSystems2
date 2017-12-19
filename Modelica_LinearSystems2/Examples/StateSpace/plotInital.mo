within Modelica_LinearSystems2.Examples.StateSpace;
function plotInital "Example plotting initial condition response"
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
This example computes the initial condition response of the system
StateSpace <i>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</i> to the initial condition <i>x0=[1;1]</i>.
</p>
</html>"));
end plotInital;
