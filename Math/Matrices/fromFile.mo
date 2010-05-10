within Modelica_LinearSystems2.Math.Matrices;
function fromFile "Read matrix from a matlab file"
  input String fileName=DataDir + "m.mat" 
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="matrix file")));
  input String matrixName="m" "Name of the matrix";

protected
  input Integer Msizes[2]=readMatrixSize(fileName, matrixName);
  input Integer n=Msizes[1];
  input Integer m=Msizes[2];
  Real M[n,m]=readMatrix(
      fileName,
      matrixName,
      n,
      m);

public
  output Real A[n,m]=M;
algorithm

end fromFile;
