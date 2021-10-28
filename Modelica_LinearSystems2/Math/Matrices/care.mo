within Modelica_LinearSystems2.Math.Matrices;
function care "Solution of continuous-time algebraic Riccati equations"
  import Modelica_LinearSystems2.Math.Matrices;
  import Complex;

  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),:];
  input Real R[size(B, 2),size(B, 2)]=identity(size(B, 2));
  input Real Q[size(A, 1),size(A, 1)]=identity(size(A, 1));
  input Boolean refine = false;
protected
  Integer n=size(A, 1);
  Real G[size(A, 1),size(A, 1)]=B*Modelica.Math.Matrices.solve2(R, transpose(B));
  Real H[2*size(A, 1),2*size(A, 1)]=[A,-G; -Q,-transpose(A)];
  Real H_RSF[2*size(A, 1),2*size(A, 1)]=H;
  Real Z[size(H, 1),size(H, 2)];
  Real Z11[size(A, 1),size(A, 2)];
  Real Z21[size(A, 1),size(A, 2)];
  Real alphaReal[size(H, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
  Real alphaImag[size(H, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";
  Integer info;
  Integer evSize;
  Complex xc[2];
public
  output Real X[size(A, 1),size(A, 2)] "Stabilizing solution of CARE";
  output Complex ev[size(A, 1)] "Eigenvalues of the closed loop system";

algorithm
  if n > 1 then
    (H_RSF,Z,alphaReal,alphaImag) := Modelica.Math.Matrices.realSchur(H);
    (H_RSF,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSF(
      true,
      H_RSF,
      Z,
      alphaReal,
      alphaImag);
    evSize := size(ev, 1);
    for i in 1:evSize loop
      ev[i] := Complex(alphaReal[i], alphaImag[i]);
    end for;

    Z11 := Z[1:n, 1:n];
    Z21 := Z[n + 1:2*n, 1:n];
    if size(Z11, 1) > 0 then
//  X := transpose(Matrices.solve2(transpose(Z11), transpose(Z21)));
      (X,info) := Matrices.LAPACK.dgesvx(Z11, transpose(Z21));
      //this function does not need to transpose Z11 as solve2 does
      X := transpose(X);
      assert(info == 0,
        "Solving a linear system of equations with function \"Matrices.LAPACK.dgesvx\" is not possible, because the system has either no or infinitely many solutions (input A is singular).");
      if refine then
        X := Matrices.Internal.carenls(A, B, R, Q, X);
      end if;
    else
      X := fill(0,size(Z21, 1),size(Z11, 1));
    end if;

  elseif n == 1 then
//    xc := Polynomial.roots(Polynomial({-G[1, 1],2*A[1, 1],Q[1, 1]}));
//    X := matrix(-abs(xc[1].re));

    X := matrix((A[1,1]-sqrt(A[1,1]*A[1,1]+G[1,1]*Q[1,1]))/G[1,1]);
    if X[1,1]*G[1,1]<A[1,1] then
      X:=matrix((A[1, 1] + sqrt(A[1, 1]*A[1, 1] + G[1, 1]*Q[1, 1]))/G[1, 1]);
    end if;
  else
    X := fill(0, 0, 0);
  end if;

  annotation (Documentation(info="<html>
<p>
Function <strong>care</strong> computes the solution <strong>X</strong> of the continuous-time
algebraic Riccati equation
</p>
<blockquote><pre>
<strong>Q</strong> + <strong>A</strong>'*<strong>X</strong> + <strong>X</strong>*<strong>A</strong> - <strong>X</strong>*<strong>G</strong>*<strong>X</strong> = <strong>0</strong>
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
<strong>G</strong> = <strong>B</strong>*<strong>R</strong><sup>-1</sup>*<strong>B</strong>'
</pre></blockquote>
<p>
using the Schur vector approach proposed by Laub [1].
</p>
<p>
It is assumed that <strong>Q</strong> is symmetric and positve semidefinite and <strong>R</strong> is
symmetric, nonsingular and positive definite, (<strong>A</strong>,<strong>B</strong>) is stabilizable
and (<strong>A</strong>,<strong>Q</strong>) is detectable.
<strong>
The assumptions are not checked in this function!
</strong>
</p>
<p>
The assumptions guarantee that Hamiltonian matrix
</p>
<blockquote><pre>
<strong>H</strong> = [<strong>A</strong>, -<strong>G</strong>; -<strong>Q</strong>, -<strong>A</strong>']
</pre></blockquote>
<p>
has no pure imaginary eigenvalue and can be put to an ordered real Schur form
</p>
<blockquote><pre>
<strong>U</strong>'*<strong>H</strong>*<strong>U</strong> = <strong>S</strong> = [<strong>S</strong>11, <strong>S</strong>12; <strong>0</strong>, <strong>S</strong>22]
</pre></blockquote>
<p>
with orthogonal similarity transformation <strong>U</strong>. <strong>S</strong> is ordered in such a way,
that <strong>S11</strong> contains the n stable eigenvalues of the closed loop system with system matrix
</p>
<blockquote><pre>
<strong>A</strong> - <strong>B</strong>*<strong>R</strong><sup>-1</sup>*<strong>B</strong>'*<strong>X</strong>
</pre></blockquote>
<p>
If <strong>U</strong> is partitioned to
</p>
<blockquote><pre>
<strong>U</strong> = [<strong>U</strong>11, <strong>U</strong>12; <strong>U</strong>21, <strong>U</strong>22]
</pre></blockquote>
<p>
with dimenstions according to <strong>S</strong>, the solution <strong>X</strong> can be calculated by
</p>
<blockquote><pre>
<strong>X</strong>*<strong>U</strong>11 = <strong>U</strong>21.
</pre></blockquote>
<p>
The algorithm uses LAPACK routines dgehrd (to compute the upper Hessenberg matrix of <strong>H</strong>), dorghr (to calculate the orthogonal
matrix from the elementary reflectors as returned from dgehrd), dhseqr (to put transformed <strong>H</strong> to Schur form and to calculate the eigenvalues
of the closed loop system) and dtrsen (to compute the ordered real Schur form and matrix <strong>U</strong>).
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Laub A.J. (1979):</dt>
<dd> <strong>A Schur Method for Solving Algebraic Riccati equations</strong>.
     IEEE Trans. Auto. Contr., AC-24, pp. 913-921.<br>&nbsp;</dd>
</dl>
</html>",
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
end care;
