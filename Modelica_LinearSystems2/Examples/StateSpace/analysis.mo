within Modelica_LinearSystems2.Examples.StateSpace;
function analysis "Example to check controllability of a state space system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;

  input StateSpace ssi = Modelica_LinearSystems2.StateSpace(
    A=[-3,2,-3,4,5,6; 0,6,7,8,9,4; 0,2,3,0,78,6; 0,1,2,2,3,3; 0,13,34,0,0,1; 0,
       0,0,-17,0,0],
    B=[1,0; 0,1; 1,0; 0,1; 1,0; 0,1],
    C=[0,0,1,0,1,0; 0,1,0,0,1,1],
    D=[0,0; 0,0],
    xNames={"x1","x2","x3","x4","x5","x6"},
    uNames={"u1","u2"}, yNames={"y1","y2"});
  input Internal.AnalyseOptions analyseOptions=
    Modelica_LinearSystems2.Internal.AnalyseOptions(
      plotEigenValues=true,
      plotInvariantZeros=true,
      plotStepResponse=true,
      plotFrequencyResponse=true,
      printEigenValues=true,
      printEigenValueProperties=true,
      printInvariantZeros=true,
      printControllability=true,
      printObservability=true,
      headingEigenValues="Eigenvalues",
      headingInvariantzeros="Invariant zeros",
      headingStepResponse="Step response",
      headingFrequencyResponse="Frequency response");

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

algorithm
  ok := false;
  StateSpace.Analysis.analysis(ss, fileName="analysis.html", analyseOptions=analyseOptions, description="Description of the system");
  ok := true;

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This example shows the usage of function
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Design.assignPolesMI\">StateSpace.Design.assignPolesMI</a>
which is to design pole assignment controllers for state space systems with multiple input.
</p>
</html>"));
end analysis;
