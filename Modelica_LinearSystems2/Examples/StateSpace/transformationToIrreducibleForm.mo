within Modelica_LinearSystems2.Examples.StateSpace;
function transformationToIrreducibleForm
  "Example to compute the minimal state space realization of a given SISO state space realization"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;

  input String fileName=DataDir + "abcd_siso2.mat"
    "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix" annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=fill(
      0,
      0,
      0);
  input Real B[:,:]=fill(
      0,
      0,
      1);
  input Real C[:,:]=fill(
      0,
      1,
      0);
  input Real D[:,:]=fill(
      0,
      1,
      1);

  output Boolean ok;

protected
  Boolean systemOnFile=fileName <> "NoName";
  StateSpace ss=if systemOnFile then
      Modelica_LinearSystems2.StateSpace.Import.fromFile(  fileName) else
      Modelica_LinearSystems2.StateSpace(
      A=A,
      B=B,
      C=C,
      D=D);
  StateSpace sso=Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(
                                                                 ss);

algorithm
  Modelica.Utilities.Streams.print("original system:\n" + String(ss));
  Modelica.Utilities.Streams.print("\n\nminimal system:\n" + String(sso));
  ok := true;
  annotation (Documentation(info="<html>
<p>
This example shows the usage of <b>function Modelica_LinearSystems2.StateSpace.reduceSystem</b> which compute
a controllable and observable state space realization of a given state space realization.
</p>
</html>"));
end transformationToIrreducibleForm;
