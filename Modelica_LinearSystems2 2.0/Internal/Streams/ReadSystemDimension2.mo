within Modelica_LinearSystems2.Internal.Streams;
function ReadSystemDimension2
  "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Internal.Streams;
  input String fileName=DataDir + "ss_siso.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  output Integer xuy[3];

protected
  Real sizeA[1,1]=Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
      fileName,
      "nx",
      1,
      1);

  Integer ABCDsizes[2]=
      Modelica_LinearSystems2.Internal.Streams.readMatrixOnFileSize(fileName,
      matrixName);

algorithm
  xuy[1] := integer(sizeA[1, 1]);
  xuy[2] := ABCDsizes[2] - xuy[1];
  xuy[3] := ABCDsizes[1] - xuy[1];

end ReadSystemDimension2;
