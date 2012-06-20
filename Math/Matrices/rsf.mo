within Modelica_LinearSystems2.Math.Matrices;
function rsf "Computes the real Schur form (RSF) of a square matrix"
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices.Internal;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real A[:,size(A, 1)];

protected
  Integer n=size(A, 1);
  Integer i;
  Integer lwork;
  Real tau[max(0, size(A, 1) - 1)];

  Real Aout[size(A, 1),size(A, 2)];
  Real H[size(A, 1),size(A, 2)] "upper Hessenberg form of A";
  Real Q[size(A, 1),size(A, 2)]
    "represents the Hessenberg transformation as a product of the elementary reflectors";
  Integer info1;
  Integer info2;

public
  output Real T[size(A, 1),size(A, 2)];
  output Real Z[size(A, 1),size(A, 2)];
  output Real alphaReal[size(A, 1)]
    "Real part of eigenvalue=alphaReal+i*alphaImag";
  output Real alphaImag[size(A, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

algorithm
  if size(A, 1) > 1 then

// dgehrd reduces a real matrix A to upper Hessenberg form Aout by an orthogonal
// similarity transformation:  Q' * A * Q = Aout. Q can be computed from
// Aout and tau (dorghr)

    (Aout,tau,info1) := LAPACK.dgehrd(
      A,
      1,
      n);
    assert(info1 == 0, "The " + String(-info1) +
      "'th argument of LAPACK.dgehrd had an illegal value");
// dorghr to compute Q

    (Q,info1) := LAPACK.dorghr(
      Aout,
      1,
      n,
      tau);
    assert(info1 == 0, "The " + String(-info1) +
      "'th argument of LAPACK.dorghr had an illegal value");

    H[1:2, :] := Aout[1:2, :];
    for i in 3:n loop
      H[i, i - 1:n] := Aout[i, i - 1:n];
    end for;
// dhseqr computes the eigenvalues of a real upper Hessenberg matrix H,
// the Schur form T of H that is also the Schur form of A as well as the matrix
// Z containing the Schur vectors to get A = Q*H*Q' = (Z)*T*(Z)'

    lwork := max(Internal.dhseqr_workdim(H), 1);
    (alphaReal,alphaImag,info2,T,Z) := LAPACK.dhseqr(
      H,
      lwork,
      false,
      "V",
      Q);
    assert(info2 == 0, "The output info of LAPACK.dhseqr should be zero, else if\n
     info < 0:  if info = -i, the i-th argument of LAPACK.dhseqr had an illegal value\n
     info > 0:  if INFO = i, LAPACK.dhseqr failed to compute all of the  
                 eigenvalues in a total of 30*n iterations;\n  
                 elements 1:n-1 and i+1:n of WR and WI contain those  
                 eigenvalues which have been successfully computed.\n");
  else
    T := A;
    if size(A, 1) > 0 then
      Z := [1];
      alphaReal := {1};
      alphaImag := {0};
    else
      Z := fill(
        1,
        0,
        0);
      alphaReal := fill(1, 0);
      alphaImag := fill(0, 0);
    end if;
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  (T, Z, alphaReal, alphaImag) </td><td align=center> =  </td>  <td> Modelica_LinearSystems2.Math.Matrices.<b<rsf</b>(A)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Function <b>rsf</b> (real Schur form) calculates the real Schur form af a real square matrix <b>A</b>, i.e.
<blockquote><pre>
         T
A = Z*T*Z

</pre></blockquote>
with the real nxn matrices <b>T</b> and <b>Z</b>. <b>Z</b> is an orthogonal matrix.  <b>T</b> is an block upper triangular matrix with 1x1 and 2x2 blocks in the diagonal.The 1x1 blocks contains the real eigenvalues of a. The 2x2 blocks are matrices with the conjugated complex pairs of eigenvalues, whereas the real parts of the eigenvalues are the elements of the diagonal.
<p>
The calculation is performed stepwise using several lapack routines. First, lapack.dgehrd reduces matrix <b>A</b> is to upper Hessenberg form <b>H</b>=<b>Q'AQ</b>, whereas <b>Q</b> is computed by lapack.dodrghr.Finally, lapack.dhseqr transforms <b>H</b> to <b>T</b>. The eigenvalues of <b>A</b> are calculated straightforward from <b>T</b>.
<p>
Function <b>rsf</b> does not apply lapack.dgees, a routine to directly compute the real Schur from. See also
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.rsf2\">Math.Matrices.rsf2</a>
</p>


<h4>Example</h4>
<blockquote><pre>
   Real A[3,3] = [1, 2, 3; 4, 5, 6; 7, 8, 9];
   Real T[3,3];
   Real Z[3,3];
   Real alphaReal[3];
   Real alphaImag[3];

<b>algorithm</b>
  (T, Z, alphaReal, alphaImag):=Modelica_LinearSystems2.Math.Matrices.rsf(A);
//   T = [16.12, 4.9,   1.59E-015;
//        0,    -1.12, -1.12E-015;
//        0,     0,    -1.30E-015]
//   Z = [-0.23,  -0.88,   0.41;
//        -0.52,  -0.24,  -0.82;
//        -0.82,   0.4,    0.41]
//alphaReal = {16.12, -1.12, -1.32E-015}
//alphaImag = {0, 0, 0}

</pre></blockquote>
</html> "));
end rsf;
