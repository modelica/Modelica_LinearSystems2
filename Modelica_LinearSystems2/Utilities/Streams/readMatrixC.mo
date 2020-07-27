within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixC "Read the output matrix of a state space system"
  extends Internal.partialReadStateSpaceMatrix;
  extends Modelica.Icons.Function;

public
  output Real C[ny,nx] = matrixABCD[nx + 1:nx + ny, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
C = Streams.<strong>ReadMatrixC</strong>(fileName, matrixName);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and reads the matrix&nbsp;C of
a&nbsp;state space system from this file.
</p>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read C matrix of the linearized system
ReadMatrixC(dslin.mat, \"ABCD\")
//  = 
// [0.9999999999670585, 0.0, 0.0, 0.0, 0.0, 0.0;
// 0.0, 1.0000000001008236, 0.0, 0.0, 0.0, 0.0;
// 0.0, 0.0, 1.0000000005428915, 0.0, 0.0, 0.0;
// 0.0, 0.0, 0.0, 1.000000000091112, 0.0, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.9999999999578305, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.0, 1.000000000038966;
// 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
// 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]</pre></blockquote>
</html>"));
end readMatrixC;
