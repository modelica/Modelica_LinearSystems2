within Modelica_LinearSystems2.Examples.StateSpace;
function analysisDcGain "Example to compute dcGain of a state space system"
  extends Modelica.Icons.Function;

  import Modelica.Math.Matrices;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica.Utilities.Streams.print;

  input StateSpace ssi = StateSpace(
    A=[1,0,0,0,0,0; 1,4,0,2,0,-1; 0,2,3,0,78,6; 1,1,2,2,3,3; 10,13,34,0,0,1; 3,
       0,0,2,0,0],
    B=[0,0; 0,0; 0,0; 0,0; 1,0; 0,0],
    C=[1,0,1,0,1,0; 0,1,0,1,0,1; 0,1,0,1,0,1],
    D=[0,0; 0,0; 0,0]);

  output Boolean ok;
protected
  Boolean finite;

  Real K1[size(ssi.C,1), size(ssi.B,2)];

  Real K2[2,1];
  StateSpace ss2 = StateSpace(
    A=[0,0;
       0,1],
    B=[0;-2],
    C=[1,0;0,1],
    D=[0;0]);

  Real K3[2,1];
  StateSpace ss3 = StateSpace(
    A=[0,0;
       0,1],
    B=[5;-2],
    C=[1,0;0,1],
    D=[0;0]);

algorithm
  ok := false;
  (K1, finite) :=StateSpace.Analysis.dcGain(ssi);
  print(Matrices.toString(K1, name="K1") + "\nfinite1 = " + String(finite));

  (K2, finite) :=StateSpace.Analysis.dcGain(ss2);
  print(Matrices.toString(K2, name="K2") + "\nfinite2 = " + String(finite) + " \n");

  (K3, finite) :=StateSpace.Analysis.dcGain(ss3);
  print(Matrices.toString(K3, name="K3") + "\nfinite3 = " + String(finite) + "\n");

  ok := true;

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This example shows how to calculate the dc-gain of a state space system.
</p>
</html>"));
end analysisDcGain;
