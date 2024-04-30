within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixA2 "Read the state matrix of a state space system"
  extends Modelica.Icons.Function;

  input String fileName=DataDir + "abcd.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";

protected
  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      nx + nu);

public
  output Real A[nx,nx]=ABCD[1:nx, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
A = Streams.<strong>readMatrixA2</strong>(fileName, matrixName, nx);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and returns the submatrix&nbsp;A of the matrix
<code>matrixName</code> given as
</p>
<blockquote><pre>
A[nx, nx] = matrixName[1:nx, 1:nx];
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read A matrix of the linearized system
readMatrixA2(\"dslin.mat\", \"ABCD\", 6)
//  = 
// [0.0, 1.0000000001008236, 0.0, 0.0, 0.0, 0.0;
// 0.0, 0.0, -2.6206364377231033, 1.2573088155408514, -0.779772901933996, 0.1247427420611989;
// 0.0, 0.0, 0.0, 1.000000000091112, 0.0, 0.0;
// 0.0, 0.0, -0.19138983233984724, 3.539452887205103, -3.1119345375319445, 1.6015423818240417;
// 0.0, 0.0, 0.0, 0.0, 0.0, 1.000000000038966;
// 0.0, 0.0, -7.229008278422499, 1.3904525205542553, 3.4263760991964918, -1.0645822524481547]
</pre></blockquote>
</html>"));
end readMatrixA2;
