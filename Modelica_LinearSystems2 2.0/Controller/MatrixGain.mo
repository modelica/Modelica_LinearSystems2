within Modelica_LinearSystems2.Controller;
block MatrixGain
  "Output the product of a gain matrix with the input signal vector. The matrix can be loaded from a file optionally"

  extends Modelica.Blocks.Interfaces.MIMO(final nin=size(K2, 2), final nout=
        size(K2, 1));
  extends Interfaces.PartialBlockIcon;

  parameter Boolean matrixOnFile=false
    "true if matrix should be read from file";
  parameter String fileName=Modelica_LinearSystems2.DataDir + "k.mat"
                              annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="matrix data file"),enable = matrixOnFile));
  parameter String matrixName="K" "Name of the matrix" annotation(Dialog(enable = matrixOnFile));

  parameter Real K[:,:]=[1] "Matrix  gain"          annotation(Dialog(enable = not matrixOnFile));

protected
  parameter Integer mn[2]=if matrixOnFile then readMatrixSize(fileName, matrixName) else size(K);
  parameter Integer m=mn[1];
  parameter Integer n=mn[2];
  parameter Real K2[:,:]=if matrixOnFile then
      Modelica_LinearSystems2.Math.Matrices.Internal.readMatrixGain(
      fileName,
      matrixName,
      m,
      n) else K;

  annotation (
    Documentation(info="<html>
<p>
This block is similar to Modelica.Blocks.Math.MatrixGain. Additionally
this block offers to load the matrix from a MATLAB-file. It
computes output vector <b>y</b> as <i>product</i> of the
gain matrix <b>K</b> with the input signal vector <b>u</b>:
</p>
<pre>
    <b>y</b> = <b>K</b> * <b>u</b>;
</pre>
<p>
Example:
</p>
<pre>
   parameter: <b>K</b> = [0.12 2; 3 1.5]

   results in the following equations:

     | y[1] |     | 0.12  2.00 |   | u[1] |
     |      |  =  |            | * |      |
     | y[2] |     | 3.00  1.50 |   | u[2] |
</pre>



</html>
"), Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Text(
          extent={{-88,-62},{92,58}},
          lineColor={160,160,164},
          textString="*[kij]"), Text(
          extent={{-150,150},{150,110}},
          textString="%name",
          lineColor={0,0,255})}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics));
equation
  y = K2*u;
end MatrixGain;
