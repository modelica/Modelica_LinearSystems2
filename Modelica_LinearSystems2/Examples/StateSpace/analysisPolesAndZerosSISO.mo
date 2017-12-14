within Modelica_LinearSystems2.Examples.StateSpace;
function analysisPolesAndZerosSISO
  "Compute poles and invariant zeros of a SISO state space system by transformation to a minmal system"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Complex;

  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
    annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"
    annotation(Dialog(group="system data definition",enable = systemOnFile));

  input Real A[:,:]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real B[:,:]=[1.0; 1.0; 0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real C[:,:]=[1.0,1.0,1.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));
  input Real D[:,:]=[0.0] annotation(Dialog(group="system matrices",enable = not systemOnFile));

  output Boolean ok;

protected
  Boolean systemOnFile=fileName <> "NoName";

  StateSpace ss=if systemOnFile then
    StateSpace.Import.fromFile(  fileName, matrixName) else
    StateSpace(A=A, B=B, C=C, D=D);
  StateSpace ssm=StateSpace.Transformation.toIrreducibleForm(ss);
  Complex poles[:];
  Complex zeros[:];

algorithm
  poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);

  Modelica.Utilities.Streams.print("\nThe poles of the unreduced system are\n");
  for i in 1:size(poles, 1) loop
    Modelica.Utilities.Streams.print("pole_" + String(i) + " = " + String(poles[i]));
  end for;

  zeros := StateSpace.Analysis.invariantZeros(ss);

  Modelica.Utilities.Streams.print("\n\nThe zeros of the unreduced system are\n");
  for i in 1:size(zeros, 1) loop
     Modelica.Utilities.Streams.print("zero_" + String(i) + " = " + String(zeros[i]));
  end for;

  poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ssm.A);

  Modelica.Utilities.Streams.print("\n\nThe poles of the reduced system are\n");
  for i in 1:size(poles, 1) loop
    Modelica.Utilities.Streams.print("pole_" + String(i) + " = " + String(poles[i]));
  end for;

  zeros := StateSpace.Analysis.invariantZeros(ssm);

  Modelica.Utilities.Streams.print("\n\nThe zeros of the reduced system are\n");
  for i in 1:size(zeros, 1) loop
     Modelica.Utilities.Streams.print("zero_" + String(i) + " = " + String(zeros[i]));
  end for;

  ok := true;

  annotation (Documentation(info="<html>
<p>
This example shows the computation of the poles and zeros of a SISO system and of its irreducible representation.
</p>
</html>"));
end analysisPolesAndZerosSISO;
