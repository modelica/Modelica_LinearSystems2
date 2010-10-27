within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixB2 "Read the input matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer nu "number of inputs";

protected
  Integer ABCDsizes[2]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixOnFileSize(fileName,
      matrixName);

  Real ABCD[ABCDsizes[1],nx + nu]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
      fileName,
      matrixName,
      ABCDsizes[1],
      nx + nu);

public
  output Real B[nx,nu]=ABCD[1:nx, nx + 1:nx + nu];
algorithm

end ReadMatrixB2;
