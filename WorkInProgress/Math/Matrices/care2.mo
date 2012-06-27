within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function care2
  "Solution of continuous-time algebraic Riccati equations with diagonal matrix R"
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Complex;

  input Real A[:,size(A, 1)];
  input Real B[size(A, 1),:];
  input Real R[size(B, 2),size(B, 2)]=identity(size(B, 2));
  input Real Q[size(A, 1),size(A, 1)]=identity(size(A, 1));
  input Boolean refine = false;
protected
  Integer n=size(A, 1);
  Real G[size(A, 1),size(A, 1)]=fill(0,size(A,1),size(A,1));//=B*Modelica.Math.Matrices.solve2(R, transpose(B));
  Real H[2*size(A, 1),2*size(A, 1)];//=[A,-G; -Q,-transpose(A)];
  Real H_RSF[2*size(A, 1),2*size(A, 1)];//=H;
  Real Z[size(H, 1),size(H, 2)];
  Real Z11[size(A, 1),size(A, 2)];
  Real Z21[size(A, 1),size(A, 2)];
  Real alphaReal[size(H, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
  Real alphaImag[size(H, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";
  Integer info;
  Integer l1;
  Integer l2;
  Integer evSize;
  Complex xc[2];
public
  output Real X[size(A, 1),size(A, 2)] "stabilizing solution of CARE";
  output Complex ev[size(A, 1)] "eigenvalues of the closed loop system";

algorithm
//create the upper triangle of symmetric matrix G
   for l1 in 1:n loop
     for l2 in l1:n loop
       for l3 in 1:size(B,2) loop
         G[l2, l1] := G[l2, l1] + B[l1, l3]*B[l2, l3]/R[l3,l3];
       end for;
     end for;
   end for;
// copy elements on and above the diagonal of G to the lower triangle
  G := symmetric(G);
  H := [A,-G; -Q,-transpose(A)];
  H_RSF := H;

  if n > 1 then
    (H_RSF,Z,alphaReal,alphaImag) := Matrices.rsf2(H);
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
     (X,info) := Matrices.solve2r(Z11, Z21);
//      (X,info) := Matrices.LAPACK.dgesvx(Z11, transpose(Z21));
//      //this function does not need to transpose Z11 as solve2 does
//      X := transpose(X);

      assert(info == 0,
        "Solving a linear system of equations with function \"Matrices.LAPACK.dgesvx\" is not possible, because the system has either no or infinitely many solutions (input A is singular).");
      if refine then
        X := Matrices.Internal.carenls(A, B, R, Q, X);
      end if;
    else
      X := fill(
        0,
        size(Z21, 1),
        size(Z11, 1));
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
 
 
Function <b>care</b> computes the solution <b>X</b> of the continuous-time algebraic Riccati equation
<blockquote><pre>
 <b>Q</b> + <b>A</b>'*<b>X</b> + <b>X</b>*<b>A</b> - <b>X</b>*<b>G</b>*<b>X</b> = <b>0</b> 
</pre></blockquote>
with
<blockquote><pre>
       -1
<b>G</b> = <b>B</b>*<b>R</b> *<b>B</b>' 
</pre>
</blockquote>
using the Schur vector approach proposed by Laub [1].
<p>
It is assumed that <b>Q</b> is symmetric and positve semidefinite and <b>R</b> is symmetric, nonsingular and positive definite,
(<b>A</b>,<b>B</b>) is stabilizable and (<b>A</b>,<b>Q</b>) is detectable.
<p><b>
The assumptions are not checked in this function
</b>
<p>
The assumptions guarantee that Hamiltonian matrix 
<blockquote><pre>
<b>H</b> = [<b>A</b>, -<b>G</b>; -<b>Q</b>, -<b>A</b>']
</pre></blockquote>
has no pure imaginary eigenvalue and can be put
to an ordered real Schur form 
<blockquote><pre>
<b>U</b>'*<b>H</b>*<b>U</b> = <b>S</b> = [<b>S</b>11, <b>S</b>12; <b>0</b>, <b>S</b>22]
</pre></blockquote>
with orthogonal similarity transformation <b>U</b>. <b>S</b> is ordered in such a way,
that <b>S11</b> contains the n stable eigenvalues of the closed loop system with system matrix
<blockquote><pre>
       -1
<b>A</b> - <b>B</b>*<b>R</b> *<b>B</b>'*<b>X</b>
</pre></blockquote>
If <b>U</b> is partitioned to 
<blockquote><pre>
<b>U</b> = [<b>U</b>11, <b>U</b>12; <b>U</b>21, <b>U</b>22]
</pre></blockquote>
with dimenstions according to <b>S</b>, the solution <b>X</b> can be calculated by
<blockquote><pre>
<b>X</b>*<b>U</b>11 = <b>U</b>21.
</pre></blockquote>
 
The algorithm uses LAPACK routines dgehrd (to compute the upper Hessenberg matrix of <b>H</b>), dorghr (to calculate the orthogonal
matrix from the elementary reflectors as returned from dgehrd), dhseqr (to put transformed <b>H</b> to Schur form and to calculate the eigenvalues 
of the closed loop system) and dtrsen (to compute the ordered real Schur form and matrix <b>U</b>).
 
<p>
 
 
<A name=\"References\"><B><FONT SIZE=\"+1\">References</FONT></B></A>
<PRE>
  [1] Laub, A.J.
      A Schur Method for Solving Algebraic Riccati equations.
      IEEE Trans. Auto. Contr., AC-24, pp. 913-921, 1979.
</PRE>
</html>",
        revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end care2;
