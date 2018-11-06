within Modelica_LinearSystems2.Examples.StateSpace;
function analysisObservability
  "Example to check controllability of a state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
    A=[1,2,3,4,5,6; 5,6,7,8,9,4; 0,2,3,0,78,6; 1,1,2,2,3,3; 10,13,34,0,0,1; 0,0,0,2,0,1],
    B=[1,0; 2,0; 0,1; 0,2; 0,0; 0,1],
    C=[1,0,1,0,1,0; 0,1,0,0,0,0],
    D=[0,0; 0,0]);

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

  Boolean isObservable;
  Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method;
  Boolean isSISO=StateSpace.Internal.isSISO(ss);
algorithm
  ok := false;

  method:=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.QR;
  isObservable := StateSpace.Analysis.isObservable(ss,method);

  if isSISO then
    if isObservable then
      Modelica.Utilities.Streams.print("pair (A, B) is observable");
    else
      Modelica.Utilities.Streams.print("pair (A, B) is not observable");
    end if;
  else
    if isObservable then
      Modelica.Utilities.Streams.print(
        "pair (A, B) is observable checked by QR factorzation");
    else
      Modelica.Utilities.Streams.print(
        "pair (A, B) is not observable checked by QR factorzation");
    end if;

    method :=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;
    isObservable := StateSpace.Analysis.isObservable(ss, method);
    if isObservable then
      Modelica.Utilities.Streams.print("pair (A, B) is observable checked by SVD");
    else
      Modelica.Utilities.Streams.print("pair (A, B) is not observable checked by SVD");
    end if;
  end if;

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows the usage of function 
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isObservable\">StateSpace.Analysis.isObservable</a>
which is to check whether a system is observable or not.
</p>
</html>"));
end analysisObservability;
