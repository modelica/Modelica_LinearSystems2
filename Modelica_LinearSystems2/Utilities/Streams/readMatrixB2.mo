within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixB2 "Read the input matrix of a state space system"
  extends Modelica.Icons.Function;

  input String fileName=DataDir + "abcd.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer nu "Number of inputs";

protected
  Integer ABCDsizes[2]=Modelica.Utilities.Streams.readMatrixSize(
    fileName, matrixName);

  Real ABCD[ABCDsizes[1],nx + nu]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      ABCDsizes[1],
      nx + nu);

public
  output Real B[nx,nu]=ABCD[1:nx, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
B = Streams.<strong>readMatrixB2</strong>(fileName, matrixName, nx, nu);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and returns the submatrix&nbsp;B of the matrix
<code>matrixName</code> given as
</p>
<blockquote><pre>
B[nx, nu] = matrixName[1:nx, nx+1:nx+nu];
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read B matrix of the linearized system
readMatrixB2(\"dslin.mat\", \"ABCD\", 6, 1)
//  = 
// [0.0;
//  0.13297862810901506;
//  0.0;
// -0.022602364424528787;
//  0.0;
// -0.11931966525935421]
</pre></blockquote>
</html>"));

end readMatrixB2;
