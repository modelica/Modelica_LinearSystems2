within Modelica_LinearSystems2.Examples.TransferFunction;
function importFromFile
  "Example showing how to read a transfer function from a file"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Polynomial;

  input Boolean systemOnFile=true
    "True, if state space system is defined on file"
    annotation(Dialog(group="System data definition"),choices(checkBox=true));

  input String fileName=DataDir + "tf_siso_1.mat"
    "File where numenator n and denominator d are stored"
    annotation(Dialog(group="System data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="transfer function data file"),enable = systemOnFile));

  input Polynomial n=Polynomial({1,2,3}) annotation(Dialog(group="Coefficients (if systemOnFile=false)",enable = not systemOnFile));
  input Polynomial d=Polynomial({4,5,6}) annotation(Dialog(group="Coefficients (if systemOnFile=false)",enable = not systemOnFile));
  output Boolean ok;

protected
  TransferFunction tf = if systemOnFile then TransferFunction.Import.fromFile(fileName) else TransferFunction(n=n, d=d);

algorithm
  Modelica.Utilities.Streams.print("TransferFunction = " + String(tf));
  ok := true;

  annotation (Documentation(info="<html>
<p>
This example shows how to import a transfer function from a Matlab file.
</p>
</html>
"));
end importFromFile;
