within Modelica_LinearSystems2.WorkInProgress.StateSpace.Transformation;
function toIrreducibleForm
  "Calculate a minimal controllable and observable block Hessenberg realization of a given SISO state-space representation "

 // test of SISO has to be added
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica;
  import Modelica_LinearSystems2.WorkInProgress;

  input StateSpace ss "State space system";

protected
  Modelica_LinearSystems2.Internal.StateSpaceR ssm1 = WorkInProgress.StateSpace.Internal.reducedCtrSystem(ss);
  Integer nx=size(ss.A, 1);
  Integer rankQ=ssm1.r;
  StateSpace ss2=StateSpace(
        A=transpose(ssm1.A[nx - rankQ + 1:nx, nx - rankQ + 1:nx]),
        B=transpose(ssm1.C[:, nx - rankQ + 1:nx]),
        C=transpose(ssm1.B[nx - rankQ + 1:nx, :]),
        D=ssm1.D);
  Integer nx2=ssm1.r;
  Modelica_LinearSystems2.Internal.StateSpaceR ssm2=
      WorkInProgress.StateSpace.Internal.reducedCtrSystem(ss2);
  Integer rankQ2=ssm2.r;
public
  output StateSpace ssm3(
    redeclare Real A[rankQ2,rankQ2],
    redeclare Real B[rankQ2,size(ss.B, 2)],
    redeclare Real C[size(ss.C, 1),rankQ2],
    redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);
algorithm
  ssm3 := StateSpace(
      A=transpose(ssm2.A[nx2 - rankQ2 + 1:nx2, nx2 - rankQ2 + 1:nx2]),
      B=transpose(ssm2.C[:, nx2 - rankQ2 + 1:nx2]),
      C=transpose(ssm2.B[nx2 - rankQ2 + 1:nx2, :]),
      D=(ssm2.D));

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toIrreducibleForm</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates a minimal controllable and observable block Hessenberg realization for a given state-space representation.
Therefore, all uncontrollable and unobservable modes are removed by performing orthogonal similarity transformations as described in [1].
</p>
<p>
This function is called to compute transfer functions of state space representations as described in [1]. Look at [1] for further details
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-4.5,  1.5,   4.0;
         -4.0,  1.0,   4.0;
         -1.5, -0.5,   1.0],
    B = [  1; 0; 1 ],
    C = [1,  0,  0],
    D = [0]);

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(ss);
//  tss=StateSpace(
      A=[-0.5],
      B=[-sqrt(0.5)],
      C=[-1/sqrt(0.5)1],
      D=[0])
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>

</html> ",
     revisions="<html>
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
end toIrreducibleForm;
