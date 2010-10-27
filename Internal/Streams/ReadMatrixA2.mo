within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixA2 "Read the state matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";

protected
  Integer ABCDsizes[2]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixOnFileSize(fileName,
      matrixName);

  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
      fileName,
      matrixName,
      nx + ny,
      nx + nu);

public
  output Real A[nx,nx]=ABCD[1:nx, 1:nx];
algorithm

end ReadMatrixA2;
