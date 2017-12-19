within Modelica_LinearSystems2.Examples.StateSpace;
function plotStep "Example plotting step response"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;

  input Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);

algorithm
 Modelica_LinearSystems2.StateSpace.Plot.step(    ss=ss);
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example computes and plots the step response of the system 
StateSpace <i>ss = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</i>.
</p>
</html>"));
end plotStep;
