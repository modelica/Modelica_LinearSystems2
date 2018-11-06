within Modelica_LinearSystems2.Examples.StateSpace;
function analysisStairCase
  "Example to check controllability of a state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
    annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"
    annotation(Dialog(group="system data definition",enable = systemOnFile));

  input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
    A=[1,2,3,4,5,6; 5,6,7,8,9,4; 0,2,3,0,78,6; 1,1,2,2,3,3; 10,13,34,0,0,1; 0,0,0,2,0,1],
    B=[1; 2; 0; 0; 0; 0],
    C=[1,0,1,0,1,0],
    D=[0]);

  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
  Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else ssi;
  StateSpace ss2 = StateSpace.Internal.transposeStateSpace(ss);
  StateSpace ss3 = ss;
  Boolean isControllable;
  Boolean isObservable;
  Modelica_LinearSystems2.Internal.StateSpaceR ssR;

algorithm
  ok := false;

  (isControllable,ssR) := StateSpace.Internal.staircaseSVD(ss);
  ss3.A := ssR.A;
  ss3.B := ssR.B;
  ss3.C := ssR.C;
  ss3.D := ssR.D;
  if isControllable then
    Modelica.Utilities.Streams.print("pair (A, B) is controllable");
  else
    Modelica.Utilities.Streams.print("pair (A, B) is not controllable");
  end if;

  Modelica.Utilities.Streams.print("Transformed system is"+String(ss3));

  (isObservable,ssR) := StateSpace.Internal.staircaseSVD(ss2);
  ss3.A := ssR.A;
  ss3.B := ssR.B;
  ss3.C := ssR.C;
  ss3.D := ssR.D;
  if isObservable then
    Modelica.Utilities.Streams.print("pair (A, C) is observable");
  else
    Modelica.Utilities.Streams.print("pair (A, C) is not observable");
  end if;
  Modelica.Utilities.Streams.print("\nTransformed dual system is:\n"+String(ss3));

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows the usage of the staircase algorithm to transform a state space system in staircase form to analyze controllability.
</p>
</html>"));
end analysisStairCase;
