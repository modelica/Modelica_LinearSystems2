within Modelica_LinearSystems2.Math.Matrices.Internal;
encapsulated function readMatrixGain "Read a matrix from mat-file"

  import Modelica;

  input String fileName="matrixGain.mat" "Name of the matrix data file"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                     caption="matrix data file")));
  input String matrixName="K" "Name of the matrix" annotation(Dialog);
  input Integer m;
  input Integer n;

public
  output Real K[m,n];

algorithm
  K := Modelica.Utilities.Streams.readRealMatrix(fileName, matrixName, m, n);

  annotation (
    obsolete = "Obsolete function - use Modelica.Utilities.Streams.readRealMatrix instead");
end readMatrixGain;
