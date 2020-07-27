within Modelica_LinearSystems2.Utilities.Streams;
function readMatrixB "Read the input matrix of a state space system"
  extends Internal.partialReadStateSpaceMatrix;
  extends Modelica.Icons.Function;

public
  output Real B[nx,nu] = matrixABCD[1:nx, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
B = Streams.<strong>ReadMatrixB</strong>(fileName, matrixName);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and reads the matrix&nbsp;B of
a&nbsp;state space system from this file.
</p>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read B matrix of the linearized system
ReadMatrixB(dslin.mat, \"ABCD\")
//  = 
// [0.0;
//  0.13297862810901506;
//  0.0;
// -0.022602364424528787;
//  0.0;
// -0.11931966525935421]
</pre></blockquote>
</html>"));
end readMatrixB;
