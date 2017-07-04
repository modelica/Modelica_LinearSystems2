within Modelica_LinearSystems2.Math.Matrices.Internal;
encapsulated function readMatrixGain "OBSOLETE- - use Modelica.Utilities.Streams.readRealMatrix instead: Read a matrix from mat-file"

  import Modelica;
  import Modelica_LinearSystems2;

  input String fileName = "matrixGain.mat" "Name of the matrix data file"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="Matrix data file")));
 input String matrixName="K" "Name of the matrix" annotation(Dialog);
 input Integer m;
 input Integer n;

public
  output Real K[m,n];

algorithm
  K := Modelica.Utilities.Streams.readRealMatrix(fileName, matrixName, m, n);

  annotation (Icon(graphics={
                   Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={238,46,47},
          lineThickness=0.5,
          pattern=LinePattern.Dash)}));
end readMatrixGain;
