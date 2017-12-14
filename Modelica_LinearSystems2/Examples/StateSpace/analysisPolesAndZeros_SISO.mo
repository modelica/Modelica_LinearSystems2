within Modelica_LinearSystems2.Examples.StateSpace;
function analysisPolesAndZeros_SISO
  "Obsolete function. Use Examples.StateSpace.analysisPolesAndZerosSISO instead"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Complex;

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"   annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[:,:]=[1.0; 1.0; 0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,:]=[1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[:,:]=[0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));

  output Boolean ok;

protected
  Boolean systemOnFile=fileName <> "NoName";

  Modelica_LinearSystems2.StateSpace ss=if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile(  fileName, matrixName) else
    Modelica_LinearSystems2.StateSpace(
      A=A,
      B=B,
      C=C,
      D=D);
  StateSpace ssm=Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(
    ss);
  Complex poles[:];
  Complex zeros[:];

algorithm
  poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ssm.A);

  for i in 1:size(poles, 1) loop

    Modelica.Utilities.Streams.print("pole_" + String(i) + " = " + String(poles[i]));
  end for;

  zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ssm);

  for i in 1:size(zeros, 1) loop
     Modelica.Utilities.Streams.print("zero_" + String(i) + " = " + String(zeros[i]));
  end for;

  ok := true;


  annotation (Icon(graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={255,0,0},
          pattern=LinePattern.Dash,
          lineThickness=0.5)}));
end analysisPolesAndZeros_SISO;
