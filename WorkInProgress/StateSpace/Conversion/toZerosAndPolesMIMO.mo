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
  annotation (overloadsConstructor=true, Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> StateSpace.Conversion.<b>toZerosAndPolesMIMO</b>(ss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a matrix of ZerosAndPoles records
 <blockquote><pre>
                 product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
        zp = k*---------------------------------------------------------
                product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
of a system from state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each input-output pair.


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



<h4>References</h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V.  </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control,
vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>

</html> ", revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end toZerosAndPolesMIMO;
