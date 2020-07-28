within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixC2 "Read the output matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer ny "number of outputs";
protected
  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

  Real ABCD[nx + ny,ABCDsizes[2]]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      ABCDsizes[2]);

public
  output Real C[ny,nx]=ABCD[nx + 1:nx + ny, 1:nx];
algorithm

end readMatrixC2;
