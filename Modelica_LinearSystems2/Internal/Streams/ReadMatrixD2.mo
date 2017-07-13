within Modelica_LinearSystems2.Internal.Streams;
function ReadMatrixD2 "Read the feed forward matrix of a state space system from MATLAB MAT file"
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
  input Integer ny "Number of outputs";

protected
  Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName, matrixName, nx + ny, nx + nu);

public
  output Real D[ny,nu] = ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
D = Streams.<strong>ReadMatrixD2</strong>(fileName, matrixName, nx, nu, ny)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the given <b>feed forward matrix&nbsp;D</b> of a state space system
of the given <b>system order&nbsp;nx</b>, the given <b>number of inputs&nbsp;nu</b>
and the given <b>number of outputs&nbsp;ny</b> from this file.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadMatrixD2;
