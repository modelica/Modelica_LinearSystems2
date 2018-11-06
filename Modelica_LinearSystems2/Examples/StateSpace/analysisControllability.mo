within Modelica_LinearSystems2.Examples.StateSpace;
function analysisControllability
  "Example to check controllability of a state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Utilities.Types.StaircaseMethod;

  input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
    A=[1,0,0,0,0,0; 1,0,0,0,0,0; 0,2,3,0,78,6; 1,1,2,2,3,3; 10,13,34,0,0,1; 0,
       0,0,2,0,0],
    B=[0,0; 0,0; 0,0; 0,0; 1,0; 0,0],
    C=[1,0,1,0,1,0; 0,1,0,1,0,1; 0,1,0,1,0,1],
    D=[0,0; 0,0; 0,0]);

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
    annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"
    annotation(Dialog(group="system data definition",enable = systemOnFile));

  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else ssi;
  StaircaseMethod method;
  Boolean isControllable;
  Boolean isSISO=StateSpace.Internal.isSISO(ss);

algorithm
  ok := false;

  method := StaircaseMethod.QR;
  isControllable := StateSpace.Analysis.isControllable(ss, method);

  if isSISO then
    if isControllable then
      Modelica.Utilities.Streams.print("pair (A, B) is controllable");
    else
      Modelica.Utilities.Streams.print("pair (A, B) is not controllable");
    end if;
  else
    if isControllable then
      Modelica.Utilities.Streams.print(
        "pair (A, B) is controllable by QR factorzation");
    else
      Modelica.Utilities.Streams.print(
        "pair (A, B) is not controllable by QR factorzation");
    end if;

    method := StaircaseMethod.SVD;
    isControllable := StateSpace.Analysis.isControllable(ss, method);
    if isControllable then
      Modelica.Utilities.Streams.print("pair (A, B) is controllable by SVD");
    else
      Modelica.Utilities.Streams.print("pair (A, B) is not controllable by SVD");
    end if;
  end if;

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows the usage of function 
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isControllable\">StateSpace.Analysis.isControllable</a>
which is to check whether a system is controllable or not.
</p>
</html>"));
end analysisControllability;
