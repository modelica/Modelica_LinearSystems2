within Modelica_LinearSystems2.Utilities.Streams;
function readSystemDimension "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams;

  input String fileName = "stateSpace.mat"
    "File containing the matrix matrixName, e.g. A.mat, dsin.txt"
    annotation (Dialog(loadSelector(
        filter="MAT files (*.mat);; All files (*.*)",
        caption="State space system data file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix";
  output Integer xuy[3] "Order of matrixName; size of u; size of y";

protected
  Real sizeA[1, 1] = Streams.readRealMatrix(fileName, "nx", 1, 1);

  Integer ABCDsizes[2] = Streams.readMatrixSize(fileName, matrixName);

algorithm
  xuy[1] := integer(sizeA[1, 1]);
  xuy[2] := ABCDsizes[2] - xuy[1];
  xuy[3] := ABCDsizes[1] - xuy[1];

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
xuy = Streams.<strong>readSystemDimension</strong>(fileName, matrixName);
</pre></blockquote>

<h4>Description</h4>
<p>
Opens the given MATLAB MAT file and reads the dimensions of state matrix <code>matrixName</code>.
Returns the order nx of the matrix and the numbers nu and ny of inputs and outputs, ordered as
</p>
<blockquote><pre>
xuy[3] = {nx, nu, ny};
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
// Generate dslin.mat of the double pendulum example first
Modelica_LinearSystems2.Utilities.Import.linearize(
  \"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\", 1.0);

// Read dimensions of state matrix of the linearized system
readSystemDimension(\"dslin.mat\", \"ABCD\")
//   = {6, 1, 8}
</pre></blockquote>
</html>"));
end readSystemDimension;
