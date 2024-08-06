within Modelica_LinearSystems2.Math.Matrices;
function fromFile "Read matrix from a matlab file"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams;

  input String fileName=DataDir + "m.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="matrix file")));
  input String matrixName="m" "Name of the matrix";

protected
  Integer Msizes[2]=Streams.readMatrixSize(fileName, matrixName);
  Integer n=Msizes[1];
  Integer m=Msizes[2];
  Real M[n,m]=Streams.readRealMatrix(
      fileName,
      matrixName,
      n,
      m);

public
  output Real A[n,m]=M;
algorithm

end fromFile;
