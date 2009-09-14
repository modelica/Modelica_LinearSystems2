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
  output Real S[size(A, 1),size(A, 2)];
  output Real QZ[size(A, 1),size(A, 2)];
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
// the Schur form S of H that is also the Schur form of A as well as the matrix
// QZ containing the Schur vectors to get A = Q*H*Q' = (QZ)*S*(QZ)'

    lwork := max(Internal.dhseqr_workdim(H), 1);
    (alphaReal,alphaImag,info2,S,QZ) := LAPACK.dhseqr(
      H,
      lwork,
      false,
      "V",
      Q);
    assert(info2 == 0, "The output info of LAPACK.dhseqr should be zero, else if\n
     info < 0:  if info = -i, the i-th argument of Slicot.sb02md had an illegal value\n
     info > 0:  if INFO = i, LAPACK.dhseqr failed to compute all of the  
                 eigenvalues in a total of 30*n iterations;\n  
                 elements 1:n-1 and i+1:n of WR and WI contain those  
                 eigenvalues which have been successfully computed.\n");
  else
    S := A;
    if size(A, 1) > 0 then
      QZ := [1];
      alphaReal := {1};
      alphaImag := {0};
    else
      QZ := fill(
        1,
        0,
        0);
      alphaReal := fill(1, 0);
      alphaImag := fill(0, 0);
    end if;
  end if;

end rsf;
