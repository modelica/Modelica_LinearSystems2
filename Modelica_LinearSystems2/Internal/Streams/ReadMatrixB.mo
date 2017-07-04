within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixB "Read the input matrix of a state space system"
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "abcd.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";

protected
  Real sizeA[1,1] = Streams.readRealMatrix(fileName, "nx", 1, 1);
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);
  Integer nx = integer(sizeA[1, 1]);
  Integer nu = ABCDsizes[2] - nx;
  Integer ny = ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);
public
  output Real B[nx,nu] = ABCD[1:nx, nx + 1:nx + nu];
algorithm

end ReadMatrixB;
