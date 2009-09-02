within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixA "Read the state matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";

protected
  Real sizeA[1,1]=Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
      fileName,
      "nx",
      1,
      1);
  Integer ABCDsizes[2]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixOnFileSize(fileName,
      matrixName);
  Integer nx=integer(sizeA[1, 1]);

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

end ReadMatrixA;
