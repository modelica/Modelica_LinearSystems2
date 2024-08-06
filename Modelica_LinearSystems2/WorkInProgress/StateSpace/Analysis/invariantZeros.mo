within Modelica_LinearSystems2.WorkInProgress.StateSpace.Analysis;
encapsulated function invariantZeros
  "Compute invariant zeros of linear state space system"

  import Modelica;
  import MatricesMSL = Modelica.Math.Matrices;
  import Complex;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices;
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
  Real normB=max(MatricesMSL.norm(ss.B, p=1),beta_small);
  Real normA=max(MatricesMSL.norm(ss.A, p=1),beta_small);

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
      (,R,,V2) := Matrices.QR(MatricesMSL.flipLeftRight(transpose([Cr,Dr])));
      Vf := MatricesMSL.flipLeftRight(V2);
      AfBf := [Ar,Br]*Vf;
      Af := AfBf[:, 1:size(Ar, 2)];
      Bf := Vf[1:size(Ar, 1), 1:size(Ar, 2)];

      (alphaReal,alphaImag,beta,,,info) := MatricesMSL.LAPACK.dggev(Af, Bf, n);
      assert(info == 0, "Failed to compute invariant zeros with function invariantZeros(..)");

      Zeros := fill(Complex(0), size(beta, 1));
      normB:=max(MatricesMSL.norm(Bf), beta_small);
      normA:=max(MatricesMSL.norm(Af, p=1), beta_small);

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
<blockquote><pre>
zeros = StateSpace.Analysis.<strong>invariantZeros</strong>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<blockquote><pre>
der(<strong>x</strong>) = <strong>A</strong>*<strong>x</strong> + <strong>B</strong>*<strong>u</strong>
    <strong>y</strong> = <strong>C</strong>*<strong>x</strong> + <strong>D</strong>*<strong>u</strong>
</pre></blockquote>
<p>
The invariant zeros of this system are defined as the variables
s  that make the Rosenbrock matrix of the system
</p>
<blockquote><pre>
| s<strong>I-A</strong>   <strong>-B</strong> |
|           |
| <strong>C</strong>       <strong>D</strong> |
</pre></blockquote>
<p>
singular.
</p>
<p>
This function applies the algorithm described in [1] where the system (<strong>A</strong>, <strong>B</strong>, <strong>C</strong>, <strong>D</strong>) is reduced to a new system (<strong>A</strong>r, <strong>B</strong>r <strong>C</strong>r, <strong>D</strong>r) with the same zeros and with <strong>D</strong>r of full rank.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[1, 1, 1;0, 1, 1;0,0,1],
    B=[1;0;1],
    C=[0,1,1],
    D=[0]);

  Complex zeros[:];

<strong>algorithm</strong>
  zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);
// zeros = {1, 0}
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Emami-Naeini, A. and Van Dooren, P. (1982):</dt>
<dd> <strong>Computation of Zeros of Linear Multivariable Systems</strong>.
     Automatica, 18, pp. 415-430.<br>&nbsp;</dd>
</dl>
</html>",revisions="<html>
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
end invariantZeros;
