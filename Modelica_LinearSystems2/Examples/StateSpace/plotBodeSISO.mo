within Modelica_LinearSystems2.Examples.StateSpace;
function plotBodeSISO
  "Constructs a zeros-and-poles transfer function from state space representation and plots the Bode diagram with automatic determination of the frequency range to plot "
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;

  input Boolean systemOnFile = false
    "True, if state space system is defined on file"
    annotation(Dialog(group="System data definition"),choices(checkBox=true));
  input String fileName = "NoName" "file where matrix [A, B; C, D] is stored"
    annotation (
      Dialog(
        group="System data definition",
        loadSelector(filter="MAT files (*.mat);; All files (*.*)", caption="State space system data file"),
      enable = systemOnFile));

  input Real A[:,size(A, 1)]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="System matrices",enable = not systemOnFile));
  input Real B[size(A, 2),:]=[0.0,1.0; 1.0,1.0; -1.0,0.0] annotation(Dialog(group="System matrices",enable = not systemOnFile));
  input Real C[:,size(A, 1)]=[0.0,1.0,1.0; 1.0,1.0,1.0] annotation(Dialog(group="System matrices",enable = not systemOnFile));
  input Real D[size(C, 1),size(B, 2)]=[1.0,0.0; 0.0,1.0] annotation(Dialog(group="System matrices",enable = not systemOnFile));

  input Integer iu=1 "Index of input (= column of B)" annotation(Dialog(group="Input/output in case of MIMO system",enable = not systemOnFile));
  input Integer iy=1 "Index of output (= row of C)" annotation(Dialog(group="Input/output in case of MIMO system",enable = not systemOnFile));
  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else
    StateSpace(A=A, B=B, C=C, D=D);

algorithm
  assert(iu <= size(ss.B, 2) and iu > 0, "Index for input u is " + String(iu) +
    " which is not in [1, " + String(size(ss.B, 2)) + "].");
  assert(iy <= size(ss.C, 1) and iy > 0, "Index for output y is " + String(iy) +
    " which is not in [1, " + String(size(ss.C, 1)) + "].");

  Modelica_LinearSystems2.StateSpace.Plot.bodeSISO(ss, iu, iy);
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the construnction of a zeros-and-poles-transfer-function 
from a SISO state space representation and plots the Bode diagrams with automatic 
determination of the frequency range to plot.
</p>
</html>"));
end plotBodeSISO;
