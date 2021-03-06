within Modelica_LinearSystems2.Math.Matrices;
function dare "Solution of discrete-time algebraic Riccati equations"
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Complex;

  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),:];
  input Real R[size(B, 2),size(B, 2)]=identity(size(B, 2));
  input Real Q[size(A, 1),size(A, 1)]=identity(size(A, 1));
  input Boolean refine = false;

protected
  Integer n=size(A, 1);
  Real G[size(A, 1),size(A, 1)]=B*MatricesMSL.solve2(R, transpose(B));
  Real AT[:,:]=transpose(A);
  Real LU[n,n];
  Integer p[n];
  Real H[2*n,2*n];
  Real H11[n,n];
  Real H12[n,n];
  Real H21[n,n];
  Real H22[n,n];
//  Real invAT[:,:]=transpose(MatricesMSL.inv(A));
//  Real H[:,:]=[A + G*invAT*Q,-G*invAT; -invAT*Q,invAT];
  Real H_RSF[2*n,2*n];
  Real Z[size(H, 1),size(H, 2)];
  Real Z11[size(A, 1),size(A, 2)];
  Real Z21[size(A, 1),size(A, 2)];
  Real alphaReal[size(H, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
  Real alphaImag[size(H, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

  Integer info;
  Integer evSize;
  Complex j = Modelica_LinearSystems2.Math.Complex.j();

public
  output Real X[size(A, 1),size(A, 2)]
    "Orthogonal matrix of the Schur vectors associated to ordered rsf";
  output Complex ev[size(A, 1)] "Eigenvalues of the closed loop system";

algorithm
  (LU,p) := MatricesMSL.LU(AT);
  H21 := MatricesMSL.LU_solve2(LU,p,-Q);
  H22 := MatricesMSL.LU_solve2(LU,p,identity(n));
  (LU,p) := MatricesMSL.LU(A);
  H12 := MatricesMSL.LU_solve2(LU,p,-G);
  H12 := transpose(H12);
  H11 := A - H12*Q;
  H := [H11, H12; H21, H22];
  (H_RSF,Z,alphaReal,alphaImag) := Matrices.rsf(H);// put H to Schur form
  (H_RSF,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSF(
    false,
    H_RSF,
    Z,
    alphaReal,
    alphaImag);// ordered Schur form

    evSize := size(ev, 1);
  for i in 1:evSize loop
    ev[i] := alphaReal[i] + j*alphaImag[i];
  end for;

  Z11 := Z[1:n, 1:n];
  Z21 := Z[n + 1:2*n, 1:n];
  if size(Z11, 1) > 0 then
//  X := transpose(Matrices.solve2(transpose(Z11), transpose(Z21)));
    (X,info) := Matrices.LAPACK.dgesvx(Z11, transpose(Z21));//function does not need to transpose Z11 as solve2 does
    X := transpose(X);
    assert(info == 0, "Solving a linear system of equations with function
\"Matrices.LAPACK.dgesvx\" is not possible, because the system has either
no or infinitely many solutions (input A is singular).");

    if refine then
      X := Matrices.Internal.darenls(A, B, R, Q, X);
    end if;
  else
    X := fill(0,size(Z21, 1),size(Z11, 1));
  end if;

  annotation (Documentation(info="<html>
<p>
Function <b>dare</b> computes the solution <b>X</b> of the discrete-time
algebraic Riccati equation
</p>
<blockquote><pre>
<b>X</b> = <b>A</b>'*<b>X</b>*<b>A</b> - <b>A</b>'*<b>X</b>*<b>B</b>*(<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)<sup>-1</sup>*<b>B</b>'*<b>X</b>*<b>A</b> + <b>Q</b>
</pre></blockquote>
<p>
using the Schur vector approach proposed by Laub [1].
</p>
<p>
It is assumed that <b>Q</b> is symmetric and positve semidefinite and <b>R</b>
is symmetric, nonsingular and positive definite, (<b>A</b>,<b>B</b>) is stabilizable
and (<b>A</b>,<b>Q</b>) is detectable.
<b>The assumptions are not checked in this function</b>
</p>

<p>
The assumptions guarantee that the Hamiltonian matrix
</p>
<blockquote><pre>
<b>H</b> = [<b>A</b><sup>-1</sup>, -<b>A</b><sup>-1</sup>*<b>G</b>; <b>Q</b>*<b>A</b><sup>-1</sup>, <b>A</b>' + <b>Q</b>*<b>A</b><sup>-1</sup>*<b>G</b> ]
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
<b>G</b> = <b>B</b>*<b>R</b><sup>-1</sup>*<b>B</b>'
</pre></blockquote>
<p>
has no eigenvalue on the unit circle and can be put
to an ordered real Schur form
</p>
<blockquote><pre>
<b>U</b>'*<b>H</b>*<b>U</b> = <b>X</b> = [<b>S11</b>, <b>S12</b>; <b>0</b>, <b>S22</b>]
</pre></blockquote>
<p>
with orthogonal similarity transformation <b>U</b>. <b>X</b> is ordered in such a way,
that <b>S11</b> contains the n stable eigenvalues of the closed loop system with system matrix
</p>
<blockquote><pre>
<b>A</b> - <b>B</b>*(<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)<sup>-1</sup>  *<b>B</b>'*<b>X</b>*<b>A</b>
</pre></blockquote>
<p>
If <b>U</b> is partitioned to
</p>
<blockquote><pre>
<b>U</b> = [<b>U11</b>, <b>U12</b>; <b>U21</b>, <b>U22</b>]
</pre></blockquote>
<p>
according to <b>X</b>, the solution <b>X</b> can be calculated by
</p>
<blockquote><pre>
<b>X</b>*<b>U11</b> = <b>U21</b>.
</pre></blockquote>

<p>
The algorithm uses LAPACK routines dgehrd (to compute the upper Hessenberg matrix
of <b>H</b>), dorghr (to calculate the orthogonal matrix from the elementary reflectors
as returned from dgehrd), dhseqr (to put transformed <b>H</b> to Schur form and to
calculate the eigenvalues of the closed loop system) and dtrsen (to compute the ordered
real Schur form and matrix <b>U</b>).
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Laub A.J. (1979):</dt>
<dd> <b>A Schur Method for Solving Algebraic Riccati equations</b>.
     IEEE Trans. Auto. Contr., AC-24, pp. 913-921.<br>&nbsp;</dd>
</dl>
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
end dare;
