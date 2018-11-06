within Modelica_LinearSystems2.Examples.StateSpace;
function conversionToZerosAndPolesSISO
  "Example to compute a zeros and poles representation from SISO state space representation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"
    annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[:,:]=[1.0; 1.0; 0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,:]=[1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[:,:]=[0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName, matrixName) else
    Modelica_LinearSystems2.StateSpace(A=A, B=B, C=C, D=D);
  Modelica_LinearSystems2.ZerosAndPoles zp;

algorithm
  zp := StateSpace.Conversion.toZerosAndPoles(ss);
  Modelica.Utilities.Streams.print(String(zp));
  Modelica.Utilities.Streams.print("ZerosAndPoles-TransferFunction = " + String(zp));
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the conversion of a SISO zeros-and-poles system into a state space system.
</p>
</html>"));
end conversionToZerosAndPolesSISO;
