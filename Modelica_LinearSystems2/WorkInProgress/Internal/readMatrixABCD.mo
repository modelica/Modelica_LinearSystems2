within Modelica_LinearSystems2.WorkInProgress.Internal;
function readMatrixABCD "Read the ABCD matrix of a state space system from MATLAB MAT file"
  import Modelica.Utilities.Streams;

  input String fileName = "dslin.mat" "File where matrixName data is stored"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="Open MATLAB MAT file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix on file";

protected
  Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension(fileName, matrixName);
  Integer nx = xuy[1];
  Integer nu = xuy[2];
  Integer ny = xuy[3];
public
  output Real ABCD[nx + ny, nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);

algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ABCD = Streams.<strong>readMatrixABCD</strong>(fileName, matrixName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
(in format v4, v6, v7, and if HDF is supported in the Modelica tool, also v7.3),
and reads the given matrix of a state space system from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end readMatrixABCD;
