within Modelica_LinearSystems2.Examples.StateSpace;
function transformation
  "Example to demonstrate the transformation to Jordan- observabilitiy- and controllability canonical form"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.StateSpace;

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix" annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[:,:]=[1.0; 1.0; 2.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,:]=[1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[:,:]=[0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));

protected
  Boolean systemOnFile=fileName <> "NoName";
  StateSpace ss = if systemOnFile then
    StateSpace.Import.fromFile(fileName, matrixName) else StateSpace(A=A, B=B, C=C, D=D);

  StateSpace ssj = StateSpace(A=A, B=B, C=C, D=D);
  StateSpace sso = StateSpace(A=A, B=B, C=C, D=D);
  StateSpace ssc = StateSpace(A=A, B=B, C=C, D=D);

  Real ctrb[:,:]=StateSpace.Analysis.controllabilityMatrix(ss);

public
  output Boolean ok;

algorithm
  ssj := StateSpace.Transformation.toDiagonalForm(ss);
  sso := StateSpace.Transformation.toObservabilityForm(ss);
  ssc := StateSpace.Transformation.toControllabilityForm(ss);

  Modelica.Utilities.Streams.print("\nOriginal state space system:");
  Modelica.Utilities.Streams.print(String(ss, 6, "ss"));
  Modelica.Math.Matrices.toString(ctrb, "controllability matrix", 6);

  Modelica.Utilities.Streams.print("\nJordan form:");
  Modelica.Utilities.Streams.print(String(ssj, 6, "ssy"));

  Modelica.Utilities.Streams.print("\nObservability canonical form:");
  Modelica.Utilities.Streams.print(String(sso, 6, "sso"));

  Modelica.Utilities.Streams.print("\nControllability canonical form:");
  Modelica.Utilities.Streams.print(String(ssc, 6, "ssc"));

  ok := true;

  annotation (Documentation(info="<html>
<p>
Example to demonstrate the transformation of a state space representation to Jordan- observabilitiy- and controllability canonical form
</p>
</html>"));
end transformation;
