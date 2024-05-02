within Modelica_LinearSystems2.Math.Matrices;
function dare "Solution of discrete-time algebraic Riccati equations"
  extends Modelica.Icons.Function;

  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica.ComplexMath.j;
  import Modelica_LinearSystems2.Math.Matrices;

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
    (X,info) := MatricesMSL.LAPACK.dgesvx(Z11, transpose(Z21));//function does not need to transpose Z11 as solve2 does
    X := transpose(X);
    assert(info == 0, "Solving a linear system of equations with function
\"LAPACK.dgesvx\" is not possible, because the system has either
no or infinitely many solutions (input A is singular).");

    if refine then
      X := Matrices.Internal.darenls(A, B, R, Q, X);
    end if;
  else
    X := fill(0,size(Z21, 1),size(Z11, 1));
  end if;

  annotation (Documentation(info="<html>
<p>
Function <strong>dare</strong> computes the solution <strong>X</strong> of the discrete-time
algebraic Riccati equation
</p>
<blockquote><pre>
<strong>X</strong> = <strong>A</strong>'*<strong>X</strong>*<strong>A</strong> - <strong>A</strong>'*<strong>X</strong>*<strong>B</strong>*(<strong>R</strong> + <strong>B</strong>'*<strong>X</strong>*<strong>B</strong>)<sup>-1</sup>*<strong>B</strong>'*<strong>X</strong>*<strong>A</strong> + <strong>Q</strong>
</pre></blockquote>
<p>
using the Schur vector approach proposed by Laub [1].
</p>
<p>
It is assumed that <strong>Q</strong> is symmetric and positive semidefinite and <strong>R</strong>
is symmetric, nonsingular and positive definite, (<strong>A</strong>,<strong>B</strong>) is stabilizable
and (<strong>A</strong>,<strong>Q</strong>) is detectable.
<strong>The assumptions are not checked in this function</strong>
</p>

<p>
The assumptions guarantee that the Hamiltonian matrix
</p>
<blockquote><pre>
<strong>H</strong> = [<strong>A</strong><sup>-1</sup>, -<strong>A</strong><sup>-1</sup>*<strong>G</strong>; <strong>Q</strong>*<strong>A</strong><sup>-1</sup>, <strong>A</strong>' + <strong>Q</strong>*<strong>A</strong><sup>-1</sup>*<strong>G</strong> ]
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
<strong>G</strong> = <strong>B</strong>*<strong>R</strong><sup>-1</sup>*<strong>B</strong>'
</pre></blockquote>
<p>
has no eigenvalue on the unit circle and can be put
to an ordered real Schur form
</p>
<blockquote><pre>
<strong>U</strong>'*<strong>H</strong>*<strong>U</strong> = <strong>X</strong> = [<strong>S11</strong>, <strong>S12</strong>; <strong>0</strong>, <strong>S22</strong>]
</pre></blockquote>
<p>
with orthogonal similarity transformation <strong>U</strong>. <strong>X</strong> is ordered in such a way,
that <strong>S11</strong> contains the n stable eigenvalues of the closed loop system with system matrix
</p>
<blockquote><pre>
<strong>A</strong> - <strong>B</strong>*(<strong>R</strong> + <strong>B</strong>'*<strong>X</strong>*<strong>B</strong>)<sup>-1</sup>  *<strong>B</strong>'*<strong>X</strong>*<strong>A</strong>
</pre></blockquote>
<p>
If <strong>U</strong> is partitioned to
</p>
<blockquote><pre>
<strong>U</strong> = [<strong>U11</strong>, <strong>U12</strong>; <strong>U21</strong>, <strong>U22</strong>]
</pre></blockquote>
<p>
according to <strong>X</strong>, the solution <strong>X</strong> can be calculated by
</p>
<blockquote><pre>
<strong>X</strong>*<strong>U11</strong> = <strong>U21</strong>.
</pre></blockquote>

<p>
The algorithm uses LAPACK routines dgehrd (to compute the upper Hessenberg matrix
of <strong>H</strong>), dorghr (to calculate the orthogonal matrix from the elementary reflectors
as returned from dgehrd), dhseqr (to put transformed <strong>H</strong> to Schur form and to
calculate the eigenvalues of the closed loop system) and dtrsen (to compute the ordered
real Schur form and matrix <strong>U</strong>).
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Laub A.J. (1979):</dt>
<dd> <strong>A Schur Method for Solving Algebraic Riccati equations</strong>.
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
