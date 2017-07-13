within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixC2 "Read the output matrix of a state space system from MATLAB MAT file"
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
  input Integer ny "Number of outputs";
protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

  Real ABCD[nx + ny,ABCDsizes[2]] = Streams.readRealMatrix(
    fileName, matrixName, nx + ny, ABCDsizes[2]);

public
  output Real C[ny,nx] = ABCD[nx + 1:nx + ny, 1:nx];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
C = Streams.<strong>ReadMatrixC2</strong>(fileName, matrixName, nx, ny)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the given <b>output matrix&nbsp;C</b> of a state space system
of the given <b>system order&nbsp;nx</b> and the given <b>number of outputs&nbsp;ny</b> from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadMatrixC2;
