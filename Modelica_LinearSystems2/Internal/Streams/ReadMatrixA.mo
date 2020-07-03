within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixA "Read the state matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";

protected
  Integer xuy[3] = Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension(
    fileName, matrixName);
  Integer nx = xuy[1];
  Integer nu = xuy[2];
  Integer ny = xuy[3];
  Real ABCD[nx + ny,nx + nu] = Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      nx + nu);

public
  output Real A[nx,nx]=ABCD[1:nx, 1:nx];
algorithm

end ReadMatrixA;
