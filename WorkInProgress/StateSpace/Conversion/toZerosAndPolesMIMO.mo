within Modelica_LinearSystems2.WorkInProgress.StateSpace.Conversion;
function toZerosAndPolesMIMO
  "Generate a zeros-and-poles representation from a MIMO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.WorkInProgress;

  input StateSpace ss "StateSpace object";

  output ZerosAndPoles zp[size(ss.C, 1),size(ss.B, 2)];

protected
  StateSpace ss_siso(
    redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
    redeclare Real B[size(ss.B, 1),1],
    redeclare Real C[1,size(ss.C, 2)],
    redeclare Real D[1,1]);

  Integer ny=size(ss.C, 1);
  Integer nu=size(ss.B, 2);

algorithm
  for ic in 1:ny loop
    for ib in 1:nu loop
      ss_siso := StateSpace(
          A=ss.A,
          B=matrix(ss.B[:, ib]),
          C=transpose(matrix(ss.C[ic, :])),
          D=matrix(ss.D[ic, ib]));
          zp[ic, ib] := Modelica_LinearSystems2.WorkInProgress.StateSpace.Conversion.toZerosAndPoles(ss_siso);
     end for;
  end for;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = StateSpace.Conversion.<b>toZerosAndPolesMIMO</b>(ss)
</pre> </blockquote>

<h4>Description</h4>
<p>
Computes a matrix of ZerosAndPoles records
</p>
<blockquote><pre>
          product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
zp = k * ---------------------------------------------------------
          product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
<p>
of a system from state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each input-output pair.
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

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);

// zp = [(s^2 + 5*s + 7)/( (s + 2)*(s + 3) ), 1/(s + 2);
         1/( (s + 2)*(s + 3) ), 1*(s + 1.38197)*(s + 3.61803)/( (s + 1)*(s + 2) )]
</pre></blockquote>
<p>
i.e.
</p>
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

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>

</html> ", revisions="<html>
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
end toZerosAndPolesMIMO;
