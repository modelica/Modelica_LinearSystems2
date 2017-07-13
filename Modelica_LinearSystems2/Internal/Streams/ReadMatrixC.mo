within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixC "Read the output matrix of a state space system from MATLAB MAT file"
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "abcd.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";

protected
  Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension(fileName, matrixName);
  Integer nx = xuy[1];
  Integer nu = xuy[2];
  Integer ny = xuy[3];
  Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);

public
  output Real C[ny,nx] = ABCD[nx + 1:nx + ny, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
C = Streams.<strong>ReadMatrixC</strong>(fileName, matrixName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the given <b>output matrix&nbsp;C</b> of a state space system from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadMatrixC;
