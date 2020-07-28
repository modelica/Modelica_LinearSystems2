within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixD "Read the feed forward matrix of a state space system"
  extends Internal.partialReadStateSpaceMatrix;
  extends Modelica.Icons.Function;

public
  output Real D[ny,nu] = matrixABCD[nx + 1:nx + ny, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
D = Streams.<strong>ReadMatrixD</strong>(fileName, matrixName);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and reads the matrix&nbsp;D of
a&nbsp;state space system from this file.
</p>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read D matrix of the linearized system
ReadMatrixD(dslin.mat, \"ABCD\")
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
end readMatrixD;
