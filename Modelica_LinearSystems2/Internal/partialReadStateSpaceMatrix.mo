within Modelica_LinearSystems2.Internal;
partial function partialReadStateSpaceMatrix "Read the ABCD matrix of the state space form of a system from MAT file"
  extends Modelica.Icons.Function;
  import Modelica.Utilities.Streams.readRealMatrix;

  input String fileName = "dslin.mat" "File where matrixName data is stored"
    annotation (
      Dialog(
        loadSelector(
          filter="MAT files (*.mat);; All files (*.*)",
          caption="Open MATLAB MAT file")));
  input String matrixName = "ABCD"
    "Name of the generalized state space system matrix on file";
protected
  Integer xuy[3]=Modelica_LinearSystems2.Utilities.Streams.readSystemDimension(fileName, matrixName);
  Integer nx = xuy[1];
  Integer nu = xuy[2];
  Integer ny = xuy[3];
  Real matrixABCD[nx + ny, nx + nu] = readRealMatrix(
    fileName, matrixName, nx + ny, nx + nu);

  annotation (
    Documentation(info="<html>
<p>
This <em>partial</em> function opens the given MATLAB MAT file and
reads the given matrix of a&nbsp;state space system from this file.
This function has no outputs, beining considered as
a&nbsp;&quot;partial&quot; function to be further extended.
Especially, operations on the protected <code>matrixABCD</code>
can be done in an extending function.
</p>


<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Internal.readSystemDimension\">readSystemDimension</a>
</p>
</html>"));
end partialReadStateSpaceMatrix;
