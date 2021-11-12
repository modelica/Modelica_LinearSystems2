within Modelica_LinearSystems2.WorkInProgress.StateSpace.Conversion;
function toTransferFunctionMIMO
  "Generate a transfer function of a MIMO system from state space representation"
  import Modelica_LinearSystems2;

  import Modelica;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.WorkInProgress;
  input StateSpace ss "StateSpace object";

  output TransferFunction tf[size(ss.C, 1),size(ss.B, 2)]
    "Matrix of transfer function objects";

protected
  ZerosAndPoles zp[:,:];
  parameter Integer m=size(ss.B, 2);
  parameter Integer p=size(ss.C, 1);

algorithm
  zp := Modelica_LinearSystems2.WorkInProgress.StateSpace.Conversion.toZerosAndPolesMIMO(ss);
  for i1 in 1:m loop
    for i2 in 1:p loop
      tf[i2, i1] := Modelica_LinearSystems2.WorkInProgress.ZerosAndPoles.Conversion.toTransferFunction(zp[i2, i1]);
    end for;
  end for;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> StateSpace.Conversion.<strong>toTransferFunctionMIMO</strong>(ss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a matrix of TransferFunction records
</p>
<blockquote><pre>
        n_i(s)     b0_i + b1_i*s + ... + bn_i*s^n
tf_i = -------- = --------------------------------
        d_i(s)     a0_i + a1_i*s + ... + an_i*s^n
</pre></blockquote>
<p>
with repetitive application of <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">Conversion.toTransferFunction</a>
</p>


<h4>Example</h4>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> StateSpace.Conversion.<strong>toTransferFunctionMIMO</strong>(ss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a matrix of TransferFunction records
</p>
<blockquote><pre>
      n(s)     b0 + b1*s + ... + bn*s^n
tf = ------ = --------------------------
      d(s)     a0 + a1*s + ... + an*s^n
</pre></blockquote>
<p>
with repetitive application of <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">Conversion.toTransferFunction</a>
</p>


<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-1.0, 0.0, 0.0;
          0.0,-2.0, 0.0;
          0.0, 0.0,-3.0],
      B = [0.0, 1.0;
           1.0, 1.0;
          -1.0, 0.0],
      C = [0.0, 1.0, 1.0;
           1.0, 1.0, 1.0],
      D = [1.0, 0.0;
           0.0, 1.0]);

<strong>algorithm</strong>
  zp:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);

// zp = [(s^2 + 5*s + 7)/(s^2 + 5*s + 6), 1/(s + 2);
         1/(s^2 + 5*s + 6), (1*s^2 + 5*s + 5)/(s^2 + 3*s + 2)]
</pre></blockquote>
i.e.
 <blockquote><pre>
           |                                                   |
           |    (s^2+5*s+7)                    1               |
           | -----------------               -----             |
           |  (s + 2)*(s + 3)                (s+2)             |
    tf  =  |                                                   |
           |        1             (s + 1.38197)*(s + 3.61803)  |
           | -------------       ----------------------------- |
           | (s + 2)*(s + 3)            (s + 1)*(s + 2)        |
           |                                                   |
</pre></blockquote>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end toTransferFunctionMIMO;
