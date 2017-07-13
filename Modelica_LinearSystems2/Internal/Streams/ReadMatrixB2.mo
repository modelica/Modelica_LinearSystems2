within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixB2 "Read the input matrix of a state space system from MATLAB MAT file"
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
  input Integer nu "Number of inputs";

protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);
  Real ABCD[ABCDsizes[1],nx + nu] = Streams.readRealMatrix(
    fileName, matrixName, ABCDsizes[1], nx + nu);

public
  output Real B[nx,nu] = ABCD[1:nx, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
B = Streams.<strong>ReadMatrixB2</strong>(fileName, matrixName, nx, nu)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the given <b>input matrix&nbsp;B</b> of a state space system
of the given <b>system order&nbsp;nx</b> and the given <b>number of inputs&nbsp;nu</b> from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadMatrixB2;
