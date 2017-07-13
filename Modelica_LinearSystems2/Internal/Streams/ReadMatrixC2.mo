within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixC2 "Read the output matrix of a state space system"
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "abcd.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer ny "number of outputs";
protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

  Real ABCD[nx + ny,ABCDsizes[2]] = Streams.readRealMatrix(
    fileName, matrixName, nx + ny, ABCDsizes[2]);

public
  output Real C[ny,nx] = ABCD[nx + 1:nx + ny, 1:nx];
algorithm

end ReadMatrixC2;
