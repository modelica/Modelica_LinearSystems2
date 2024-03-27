within Modelica_LinearSystems2.Utilities.Streams;
function ReadSystemDimension
  "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs"

  input String fileName=DataDir + "ss_siso.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  output Integer xuy[3];

protected
  Real sizeA[1,1]=Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      "nx",
      1,
      1);

  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

algorithm
  xuy[1] := integer(sizeA[1, 1]);
  xuy[2] := ABCDsizes[2] - xuy[1];
  xuy[3] := ABCDsizes[1] - xuy[1];

end ReadSystemDimension;
