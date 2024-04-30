within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixC2 "Read the output matrix of a state space system"
  extends Modelica.Icons.Function;

  input String fileName=DataDir + "abcd.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer ny "Number of outputs";
protected
  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

  Real ABCD[nx + ny,ABCDsizes[2]]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      ABCDsizes[2]);

public
  output Real C[ny,nx]=ABCD[nx + 1:nx + ny, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
C = Streams.<strong>readMatrixC2</strong>(fileName, matrixName, nx, ny);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and returns the submatrix&nbsp;C of the matrix
<code>matrixName</code> given as
</p>
<blockquote><pre>
C[ny, nx] = matrixName[nx+1:nx+ny, 1:nx];
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read C matrix of the linearized system
readMatrixC2(\"dslin.mat\", \"ABCD\", 6, 8)
//  = 
// [0.9999999999670585, 0.0, 0.0, 0.0, 0.0, 0.0;
// 0.0, 1.0000000001008236, 0.0, 0.0, 0.0, 0.0;
// 0.0, 0.0, 1.0000000005428915, 0.0, 0.0, 0.0;
// 0.0, 0.0, 0.0, 1.000000000091112, 0.0, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.9999999999578305, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.0, 1.000000000038966;
// 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
</pre></blockquote>
</html>"));
end readMatrixC2;
