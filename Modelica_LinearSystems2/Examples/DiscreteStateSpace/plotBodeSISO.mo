within Modelica_LinearSystems2.Examples.DiscreteStateSpace;
function plotBodeSISO "Example to plot Bode diagram of a discrete state space system"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  input DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A=[0.995166584721977,0.0950040833529266,0.0,0.0,0.0,0.0; -0.0950040833529266,
       0.90016250136905,0.0,0.0,0.0,0.0; 0.00676921321653657,0.00700084602310958,
       0.98598654522888,0.0901824307998049,0.0,0.0; 0.128272800176598,0.135042013393134,
       -0.270547292399415,0.80562168362927,0.0,0.0; 0.00349754076724123,0.0036162394148413,
       0.0879241048980837,0.0473495413016237,0.90483741803596,0.0; 0.000171498804565003,
       0.000175864915788407,0.00655179215601711,0.00344328523986849,0.132388324832705,
       0.860707976425058],
    B=[0.00483341527802304; 0.0950040833529266; 0.0072442415545839; 0.142274492222817;
       0.00374093629871555; 0.000180407781654784],
    C=[0.0,0.0,0.0,0.0,0.0,0.888888888888889],
    D=[0.0],
    Ts=0.1,
    B2=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0],
    method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact);
  input Integer iu=1 "index of inout";
  input Integer iy=1 "index of output";
  output Boolean ok;

algorithm
  assert(iu <= size(dss.B, 2) and iu > 0, "index for input is " + String(iu) + " which is not in [1, "
     + String(size(dss.B, 2)) + "].");
  assert(iy <= size(dss.C, 1) and iy > 0, "index for output is " + String(iy) + " which is not in [1, "
     + String(size(dss.C, 1)) + "].");

  Modelica_LinearSystems2.DiscreteStateSpace.Plot.bodeSISO(dss, iu,  iy);
  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates the construnction of a zeros-and-poles transfer function from 
a SISO state space representation and plots the Bode diagrams with automatic determination 
of the frequency range to plot.
</p>
</html>"));
end plotBodeSISO;
