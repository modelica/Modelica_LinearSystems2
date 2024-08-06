within Modelica_LinearSystems2.Examples.StateSpace;
function plotRamp "Plot ramp response"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.StateSpace ss=
      Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);

algorithm
  Modelica_LinearSystems2.StateSpace.Plot.ramp(ss=ss);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the ramp response of the system
StateSpace <em>ss = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</em>.
</p>
</html>"));
end plotRamp;
