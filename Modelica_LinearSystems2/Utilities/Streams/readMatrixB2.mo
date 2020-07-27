within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixB2 "Read the input matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer nu "number of inputs";

protected
  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

  Real ABCD[ABCDsizes[1],nx + nu]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      ABCDsizes[1],
      nx + nu);

public
  output Real B[nx,nu]=ABCD[1:nx, nx + 1:nx + nu];
algorithm

end readMatrixB2;
