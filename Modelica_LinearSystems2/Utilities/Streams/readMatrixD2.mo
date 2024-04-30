within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixD2 "Read the feed forward matrix of a state space system"
  extends Modelica.Icons.Function;

  input String fileName=DataDir + "abcd.mat"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
  input String matrixName="ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";
  input Integer nu "Number of inputs";
  input Integer ny "Number of outputs";

protected
  Real ABCD[nx + ny,nx + nu]=
      Modelica.Utilities.Streams.readRealMatrix(
      fileName,
      matrixName,
      nx + ny,
      nx + nu);

public
  output Real D[ny,nu]=ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
D = Streams.<strong>readMatrixD2</strong>(fileName, matrixName, nx, nu, ny);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and returns the submatrix&nbsp;D of the matrix
<code>matrixName</code> given as
</p>
<blockquote><pre>
D[ny, nu] = ABCD[nx+1:nx+ny, nx+1:nx+nu];
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read D matrix of the linearized system
readMatrixD2(\"dslin.mat\", \"ABCD\", 6, 1, 8)
//  = 
// [0.0;
// 0.0;
// 0.0;
// 0.0;
// 0.0;
// 0.0;
// 0.0;
// 0.0]
</pre></blockquote>
</html>"));
end readMatrixD2;
