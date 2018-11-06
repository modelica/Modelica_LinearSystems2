within Modelica_LinearSystems2.Examples.StateSpace;
function plotPolesAndZeros
  "Example for plotting eigenvalues and invariant zeros of a state space system"

  import Modelica_LinearSystems2;
  extends Modelica_LinearSystems2.StateSpace.Plot.polesAndZeros(ss=
    Modelica_LinearSystems2.StateSpace(
      A=[-3,2,-3,4,5,6; 0,6,7,8,9,4; 0,2,3,0,78,6; 0,1,2,2,3,3; 0,13,34,0,0,1; 0,0,0,-17,0,0],
      B=[1,0; 0,1; 1,0; 0,1; 1,0; 0,1],
      C=[0,0,1,0,1,0; 0,1,0,0,1,1],
      D=[0,0; 0,0]));
  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the plotting of eigenvalues and invariant zeros of a state space system.
</p>
</html>"));
end plotPolesAndZeros;
