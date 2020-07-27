within Modelica_LinearSystems2.Utilities.Streams;
function readSystemDimension "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs"
  import Modelica.Utilities.Streams;

  input String fileName = "stateSpace.mat"
    "File containing the matrix matrixName, e.g. A.mat, dsin.txt"
    annotation (Dialog(loadSelector(
        filter="MAT files (*.mat);; All files (*.*)",
        caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";
  output Integer xuy[3] "Order of matrixName; size of u; size of y";

protected
  Real sizeA[1, 1] = Streams.readRealMatrix(fileName, "nx", 1, 1);

  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

algorithm
  xuy[1] := integer(sizeA[1, 1]);
  xuy[2] := ABCDsizes[2] - xuy[1];
  xuy[3] := ABCDsizes[1] - xuy[1];

end readSystemDimension;
