within Modelica_LinearSystems2.Internal.Streams;
function ReadSystemDimension
  "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs from MATLAB MAT file"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica.Utilities.Streams;

  input String fileName = DataDir + "ss_siso.mat"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";
  output Integer xuy[3];

protected
  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

algorithm
  xuy[1] := integer(scalar(Streams.readRealMatrix(fileName, "nx", 1, 1)));
  xuy[2] := ABCDsizes[2] - xuy[1];
  xuy[3] := ABCDsizes[1] - xuy[1];

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
xuy = Streams.<strong>ReadSystemDimension</strong>(fileName, matrixName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function opens the given MATLAB MAT file
and reads the sizes of the given <b>matrix</b> of a state space system from this file.
Afterwards, the <b>system order&nbsp;nx</b>, the <b>number of inputs&nbsp;nu</b>
and the <b>number of outputs&nbsp;ny</b> are returned all collected in one vector <b>xuy</b>.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">Modelica.Utilities.Streams.readMatrixSize</a>,
<a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">Modelica.Utilities.Streams.readRealMatrix</a>
</p>
</html>"));
end ReadSystemDimension;
