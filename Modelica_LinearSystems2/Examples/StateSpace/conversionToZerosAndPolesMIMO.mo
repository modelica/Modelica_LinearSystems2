within Modelica_LinearSystems2.Examples.StateSpace;
function conversionToZerosAndPolesMIMO
  "Example to compute a zeros-and-poles representation of a MIMO system from state space representation"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix" annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,size(A, 1)]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[size(A, 2),:]=[0.0,1.0; 1.0,1.0; -1.0,0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,size(A, 1)]=[0.0,1.0,1.0; 1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[size(C, 1),size(B, 2)]=[1.0,0.0; 0.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));

  output Boolean ok;

protected
  Modelica_LinearSystems2.StateSpace ss = if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile( fileName, matrixName) else
    Modelica_LinearSystems2.StateSpace(A=A, B=B, C=C, D=D);

protected
  Modelica_LinearSystems2.ZerosAndPoles zp[:,:];

algorithm
  zp := Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPolesMIMO(ss);
  for i1 in 1:size(ss.C, 1) loop
    for i2 in 1:size(ss.B, 2) loop
      Modelica.Utilities.Streams.print("ZerosAndPoles-TransferFunction[" + String(i1) + "," +
        String(i2) + "] = " + String(zp[i1, i2]));
    end for;
  end for;
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the conversion of a MIMO zeros-and-poles system into a state space system.
</p>
</html>"));
end conversionToZerosAndPolesMIMO;
