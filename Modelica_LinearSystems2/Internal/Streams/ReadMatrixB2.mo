within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixB2 "Read the input matrix of a state space system"
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "abcd.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer nu "Number of inputs";

protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);
  Real ABCD[ABCDsizes[1],nx + nu] = Streams.readRealMatrix(
    fileName, matrixName, ABCDsizes[1], nx + nu);

public
  output Real B[nx,nu] = ABCD[1:nx, nx + 1:nx + nu];
algorithm

end ReadMatrixB2;
