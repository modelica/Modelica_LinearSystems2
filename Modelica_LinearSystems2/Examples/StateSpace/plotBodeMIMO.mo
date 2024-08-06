within Modelica_LinearSystems2.Examples.StateSpace;
function plotBodeMIMO
  "Constructs zeros-and-poles-transfers from state space representation and plots the Bode diagrams with automatic determination of the frequency range to plot "
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file"),enable = systemOnFile));

  input Real A[:,size(A, 1)]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0]  annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[size(A, 2),:]=[0.0,1.0; 1.0,1.0; -1.0,0.0]  annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,size(A, 1)]=[0.0,1.0,1.0; 1.0,1.0,1.0]  annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[size(C, 1),size(B, 2)]=[1.0,0.0; 0.0,1.0]  annotation(Dialog(group="system matrices",enable = not systemOnFile));

  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
    StateSpace.Import.fromFile(fileName) else
    StateSpace(A=A, B=B, C=C, D=D);

algorithm
  StateSpace.Plot.bodeMIMO(ss, onFile=true);
  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates the construnction of a zeros-and-poles-transfer-function-matrix 
from a MIMO state space representation and plots the Bode diagrams with automatic 
determination of the frequency range to plot.
</p>
</html>"));
end plotBodeMIMO;
