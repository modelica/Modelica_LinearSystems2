within Modelica_LinearSystems2.Examples.StateSpace;
function plotRamp "Example plotting ramp response"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.StateSpace ss=
      Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);

algorithm
  Modelica_LinearSystems2.StateSpace.Plot.ramp(ss=ss);

  annotation (
    Documentation(info="<html>
<p>
This example computes the ramp response of the system
StateSpace <i>ss = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</i>.
</p>
</html>"));
end plotRamp;
