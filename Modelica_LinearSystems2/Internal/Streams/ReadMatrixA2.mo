within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixA2 "Read the state matrix of a state space system from MATLAB MAT file"
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "abcd.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";
  input Integer nx "System order";

protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

  Integer nu = ABCDsizes[2] - nx;
  Integer ny = ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);

public
  output Real A[nx,nx] = ABCD[1:nx, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
A = Streams.<strong>ReadMatrixA2</strong>(fileName, matrixName, nx)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the given <b>state matrix&nbsp;A</b> of a state space system
of the given <b>system order&nbsp;nx</b> from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadMatrixA2;
