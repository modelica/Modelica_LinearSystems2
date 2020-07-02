within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixC "Read the output matrix of a state space system"
  input String fileName=DataDir + "abcd.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";

protected
  Real sizeA[1,1]=Modelica.Utilities.Streams.readRealMatrix(
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
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      nx + nu);
public
  output Real C[ny,nx]=ABCD[nx + 1:nx + ny, 1:nx];
algorithm

end ReadMatrixC;
