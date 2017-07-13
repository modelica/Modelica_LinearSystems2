within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixD2 "Read the feed forward matrix of a state space system"
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
  input Integer nu "number of inputs";
  input Integer ny "number of outputs";

protected
  Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);

public
  output Real D[ny,nu] = ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
algorithm

end ReadMatrixD2;
