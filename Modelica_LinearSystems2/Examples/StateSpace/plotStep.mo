within Modelica_LinearSystems2.Examples.StateSpace;
function plotStep "Step plot example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;

  input StateSpace ss=StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);

algorithm
  StateSpace.Plot.step(ss=ss);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes and plots the step response of the system 
StateSpace <em>ss = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</em>.
</p>
</html>"));
end plotStep;
