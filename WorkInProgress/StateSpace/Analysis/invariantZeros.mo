within Modelica_LinearSystems2.WorkInProgress.StateSpace.Analysis;
encapsulated function invariantZeros
  "Compute invariant zeros of linear state space system"

  import Modelica;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  import Modelica_LinearSystems2.WorkInProgress;

  input StateSpace ss "Linear system in state space form";

  output Complex Zeros[:]
    "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

protected
  Integer n=10;
  Integer m;
  Integer p;
  Real Ar[:,:];
  Real Br[:,:];
  Real Cr[:,:];
  Real Dr[:,:];

  Real Af[:,:];
  Real Bf[:,:];
  Real AfBf[:,:];

  Real V2[:,:];
  Real Vf[:,:];
  Real R[:,:];

  Integer na;
 Real alphaReal[:];
  Real alphaImag[:];
  Real beta[:];
  Integer info;
  Real beta_small=100*Modelica.Constants.eps;
  Real normB=max(Modelica.Math.Matrices.norm(ss.B, p=1),beta_small);
  Real normA=max(Modelica.Math.Matrices.norm(ss.A, p=1),beta_small);

algorithm
  if min(size(ss.B)) == 0 or min(size(ss.C)) == 0 then
    Zeros := fill(Complex(0), 0);
  else
    (Ar,Br,Cr,Dr,n,m,p) := WorkInProgress.StateSpace.Internal.reduceRosenbrock(ss.A, ss.B, ss.C, ss.D);
    if n > 0 then
      (Ar,Br,Cr,Dr,n,m,p) := WorkInProgress.StateSpace.Internal.reduceRosenbrock(transpose(Ar), transpose(Cr), transpose(Br), transpose(Dr));
    end if;
    if n == 0 then
      Zeros := fill(Complex(0), 0);
    else
      (,R,,V2) := Matrices.QR(Matrices.fliplr(transpose([Cr,Dr])));
      Vf := Matrices.fliplr(V2);
      AfBf := [Ar,Br]*Vf;
      Af := AfBf[:, 1:size(Ar, 2)];
      Bf := Vf[1:size(Ar, 1), 1:size(Ar, 2)];

      (alphaReal,alphaImag,beta,,,info) := LAPACK.dggev(Af, Bf, n);
      assert(info == 0, "Failed to compute invariant zeros with function invariantZeros(..)");

      Zeros := fill(Complex(0), size(beta, 1));
      normB:=max(Modelica.Math.Matrices.norm(Bf), beta_small);
      normA:=max(Modelica.Math.Matrices.norm(Af, p=1), beta_small);

// If beta[i] is zero, then zero i is infinite.
      for i in 1:size(beta, 1) loop
        if beta[i] >= normB*1e-3 then
     // finite eigenvalue
          Zeros[i].re := if abs(alphaReal[i]) >= normB*1e-12 then alphaReal[i]/
            beta[i] else 0;
          Zeros[i].im := if abs(alphaImag[i]) >= normB*1e-12 then alphaImag[i]/
            beta[i] else 0;
        end if;
      end for;

    end if;
  end if;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  zeros </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>invariantZeros</b>(ss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
s  that make the Rosenbrock matrix of the system
</p>
<pre>
    | s<b>I-A</b>   <b>-B</b> |
    |           |
    | <b>C</b>       <b>D</b> |

</pre>
singular.
<p>
This function applies the algorithm described in [1] where the system (<b>A</b>, <b>B</b>, <b>C</b>, <b>D</b>) is reduced to a new system (<b>A</b>r, <b>B</b>r <b>C</b>r, <b>D</b>r) with the same zeros and with <b>D</b>r of full rank.
</p>


</p>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[1, 1, 1;0, 1, 1;0,0,1],
      B=[1;0;1],
      C=[0,1,1],
      D=[0]);

   Complex zeros[:];

<b>algorithm</b>
  zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);
// zeros = {1, 0}

</pre></blockquote>

<h4>References</h4>
<table>
<tr> <td align=right>  [1] </td><td align=center>  Emami-Naeini, A. and Van Dooren, P. </td>  <td> \"Computation of Zeros of Linear Multivariable Systems\"  </td> <td> Automatica, 18, pp. 415-430, 1982. </td></tr>
</table>
</html> ",
         revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end invariantZeros;
