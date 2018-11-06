within Modelica_LinearSystems2.Examples.StateSpace;
function conversionToTransferFunctionSISO
  "Example to compute a transfer function from SISO state space representation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix" annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[:,:]=[1.0; 1.0; 0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,:]=[1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[:,:]=[0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  output Boolean ok;

protected
  StateSpace ss=if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile( fileName) else
    Modelica_LinearSystems2.StateSpace(A=A, B=B, C=C, D=D);
  TransferFunction tf;

algorithm
  tf := StateSpace.Conversion.toTransferFunction(ss);
  Modelica.Utilities.Streams.print("StateSpace = " + String(ss)+"\n");
  Modelica.Utilities.Streams.print("TransferFunction = " + String(tf));
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the conversion of a SISO transfer function  system into a state space system.
</p>
</html>"));
end conversionToTransferFunctionSISO;
