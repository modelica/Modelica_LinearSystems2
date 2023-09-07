within Modelica_LinearSystems2.Examples.StateSpace;
function plotBodeSISO
  "Constructs a zeros-and-poles transfer function from state space representation and plots the Bode diagram with automatic determination of the frequency range to plot "
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file";
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
    annotation (Dialog(enable=systemOnFile));

  input Real A[:,size(A, 1)]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0];
  input Real B[size(A, 2),:]=[0.0,1.0; 1.0,1.0; -1.0,0.0];
  input Real C[:,size(A, 1)]=[0.0,1.0,1.0; 1.0,1.0,1.0];
  input Real D[size(C, 1),size(B, 2)]=[1.0,0.0; 0.0,1.0];

  input Integer iu=1 "Index of inout";
  input Integer iy=1 "Index of output";
  output Boolean ok;

protected
  StateSpace ss=if systemOnFile then Modelica_LinearSystems2.StateSpace.Import.fromFile(
      fileName) else
      StateSpace(
      A=A,
      B=B,
      C=C,
      D=D);

algorithm
  assert(iu <= size(ss.B, 2) and iu > 0, "index for input is " + String(iu) +
    " which is not in [1, " + String(size(ss.B, 2)) + "].");
  assert(iy <= size(ss.C, 1) and iy > 0, "index for output is " + String(iy) +
    " which is not in [1, " + String(size(ss.C, 1)) + "].");

  Modelica_LinearSystems2.StateSpace.Plot.bodeSISO(
    ss,
    iu,
    iy);
  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates the construnction of a zeros-and-poles-transfer-function 
from a SISO state space representation and plots the Bode diagrams with automatic 
determination of the frequency range to plot.
</p>
</html>"));
end plotBodeSISO;
