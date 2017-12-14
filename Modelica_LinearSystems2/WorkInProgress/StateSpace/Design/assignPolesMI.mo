within Modelica_LinearSystems2.WorkInProgress.StateSpace.Design;
encapsulated function assignPolesMI
  "Pole assignment design algorithm for multi input systems"

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
  import Complex;

  input StateSpace ss "state space system";

  input Complex gamma[:]=fill(Complex(0), 0) "Designed Poles";
//  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
  input Real alpha=-1e10
    "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
  input Real tolerance=Modelica.Math.Matrices.norm(ss.A, 1)*1e-12
    "The tolerance to be used in determining the controllability of (A,B)";
  input Boolean calculateEigenvectors=false
    "Calculate the eigenvectors X of the closed loop system when true";

  output Real K[size(ss.B, 2),size(ss.A, 1)]
    "State feedback matrix assigning the desired poles";
  output Real S[:,:] "Closed loop System matrix";
  output Complex po[size(ss.A, 1)] "poles of the closed loop system";
  output Integer nfp
    "number of eigenvalues that are not modified with respect to alpha";
  output Integer nap "number of assigned eigenvalues";
  output Integer nup "number of uncontrollable eigenvalues";
  output Complex X[size(ss.A, 1),size(ss.A, 1)]
    "eigenvectors of the closed loop system";

protected
  Real A_rsf[size(ss.A, 1),size(ss.A, 2)];
  Real B_rsf[size(ss.B, 1),size(ss.B, 2)];
  Real Q[size(ss.A, 1),size(ss.A, 1)];
  Real Ks1[:,:];
  Real Ks2[:,:];
  Real Q2[:,:];
  Real A_rsf_1[:,:];
  Real Q1[:,:];
  Boolean select[:];
  Boolean rselectA[:];
  Real Z[:,:] "orthogonal transformation matrix";
  Real ZT[:,:] "orthogonal transformation matrix";
  Complex pf[:];
  Complex gammaReordered[:]=gamma;
  Integer info;
  Real wr[size(gamma, 1)];
  Real wi[size(gamma, 1)];
  Boolean imag=false;
  Integer i;
  Integer ii;
  Integer iii;
  Integer counter;
  Integer counter2;
  Integer n=size(ss.A, 1);
  Integer nccA "number of conjugated complex pole pairs of openloop system";
  Integer nccg "number of conjugated complex pole pairs of gamma";
  Integer rpg "number of real poles in gamma";
  Integer rpA "number of real poles of open loop system";
  Integer ncc "Min(nccA, nccg)";
  Integer rp "Min(rpg, rpA)";
  Integer ng=size(gamma,1);
  Integer nr "Differenz between rpA and rpg; Sign(rpA-rpg)*(rpA-rpg)";

  Real alphaReal[size(ss.A, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
  Real alphaImag[size(ss.A, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

  Complex SS[:,:];
  Complex Xj[:,:];
  Complex h;

  Real dist;
  Real evImag;

algorithm
  assert(size(gamma, 1) <= size(ss.A, 1),
    "At most n (order of ss) eigenvalues can be assigned");

 /* Extraction of Poles (Variable conversation) and pole sequence check */
  for i in 1:size(gamma, 1) loop
    wr[i] := gamma[i].re;
    wi[i] := gamma[i].im;
    if imag then
      assert(wi[i - 1] == -wi[i] and wr[i - 1] == wr[i],
        "Poles are in wrong sequence");
      imag := false;
    elseif wi[i] <> 0 then
      imag := true;
    end if;
  end for;

  // put matrix ss.A to real Schur form A <- QAQ' and compute B <- QB
  (A_rsf,Z,alphaReal,alphaImag) := Matrices.rsf2(ss.A);
  ZT := transpose(Z);

  // reorder real Schur form according to alpha
  (A_rsf,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSFc(
      A_rsf,
      identity(size(A_rsf, 1)),
      alphaReal,
      alphaImag,
      alpha);
  ZT := transpose(Z)*ZT;
  B_rsf := ZT*ss.B;

  // determine number of poles not to be assigned according to alpha
  nfp := 0;
  for i in 1:n loop
    if alphaReal[i] < alpha then
      nfp := nfp + 1;
    end if;
  end for;
  nap := n - nfp;

  assert(size(gamma, 1) >= nap, String(nap) +
    " poles should be modified, therefore gamma should contain at at least "
     + String(nap) + " assigned eigenvalues");

  // second reorder (reorderRSF3) according to conjugated complex pairs in A and p
  // count numbre of conjugated complex pole pairs = max(number_ccpp(eig(A), number_ccpp(gamma))
  nccA := 0;
  //mark the real poles of original system
  rselectA := fill(true, nap);
  ii := 1;
  for i in nfp + 1:n loop
    if abs(alphaImag[i]) > 0 then
      nccA := nccA + 1;
    else
      rselectA[ii] := false;

    end if;
    ii := ii + 1;
    end for;
    rpA := n-nccA;
    nccA := div(nccA, 2);

  // reorder gamma and A_rsf
  (gammaReordered,rpg) := Modelica_LinearSystems2.Internal.reorderZeros(gamma);
  gammaReordered :=Modelica.ComplexMath.Vectors.reverse(gammaReordered);
  nccg := div(size(gammaReordered, 1) - rpg, 2);
  ncc := min(nccA, nccg);
  rp := min(rpA, rpg);
  if nccA > 0 then
    (A_rsf[nfp + 1:n, nfp + 1:n],Q2) := Matrices.LAPACK.dtrsen(
        "E",
        "V",
        rselectA,
        A_rsf[nfp + 1:n, nfp + 1:n],
        identity(n - nfp));//The Schur vector matrix is identity, since A_rsf already has Schur form

    A_rsf[1:nfp, nfp + 1:n] := A_rsf[1:nfp, nfp + 1:n]*Q2;
    B_rsf[nfp + 1:n, :] := transpose(Q2)*B_rsf[nfp + 1:n, :];
    ZT[nfp + 1:n, :] := transpose(Q2)*ZT[nfp + 1:n, :];
  end if;

  // main algorithm
  K := zeros(size(ss.B, 2), size(ss.A, 1));
  counter := nfp + 1;
  counter2 := 1;

  for i in 1:rp loop // 1x1 blocks; real system pole and real assigned poles; take the next eigenvalue in the
                       // diagonal of the Schur form and search the nearest pole in the set of the real poles to assign
      dist := Modelica.Constants.inf;
      for ii in i:rpg loop // looking for nearest pole and reorder gamma
        if abs(A_rsf[n, n] - gammaReordered[ng - ii + 1].re) < dist then
          iii := ng - ii + 1;
          dist := abs(A_rsf[n, n] - gammaReordered[ng - ii + 1].re);
        end if;
      end for;
      h := gammaReordered[ng - i + 1];
      gammaReordered[ng - i + 1] := gammaReordered[iii];
      gammaReordered[iii] := h;

//###############################################################

Modelica.Utilities.Streams.print("1x1 Ann = "+String(A_rsf[n,n])+"\n ap = "+String(gammaReordered[ng - i + 1]));

//###############################################################

      Ks1 := Internal.assignOneOrTwoPoles(
        matrix(A_rsf[n, n]),
        transpose(matrix(B_rsf[n, :])),
        {gammaReordered[ng - i + 1]},
        tolerance);
      K := K + [zeros(size(Ks1, 1), size(K, 2) - 1),Ks1]*ZT;
      A_rsf := A_rsf - B_rsf*[zeros(size(Ks1, 1), size(K, 2) - 1),Ks1];
      select := fill(false, n - counter + 1);
      select[n - counter + 1] := true;

      (A_rsf[counter:n, counter:n],Q1) := Matrices.LAPACK.dtrsen(
        "E",
        "V",
        select,
        A_rsf[counter:n, counter:n],
        identity(n - counter + 1));//The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q1;
      B_rsf[counter:n, :] := transpose(Q1)*B_rsf[counter:n, :];
      ZT[counter:n, :] := transpose(Q1)*ZT[counter:n, :];
      counter := counter + 1;
      counter2 := counter2 + 1;
    end for;

    if counter2<rpg and counter2>rpA then //System has less real eigenvalues than real assigned poles
    for i in 1:div(rpg - rpA, 2) loop // 2x2 blocks; complex pair of system poles and 2 real assigned poles; take the next complex pair
                                      // (Schur bump) in the diagonal of the Schur form and search the two nearest poles in the set of the
                                      // remaining real assigned poles
      dist := Modelica.Constants.inf;
      evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);//positive imaginary part of the complex system pole pair
      for ii in 2*(i - 1) + 1:2:rpg - rpA loop
        if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag < dist then
          iii := ng - rp - ii + 1;
          dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag;
        end if;
      end for;
      h := gammaReordered[ng - rp - 2*(i - 1)];
      gammaReordered[ng - rp - 2*(i - 1)] := gammaReordered[iii];
      gammaReordered[iii] := h;
      dist := Modelica.Constants.inf;
      for ii in 2*(i - 1) + 1:2:rpg - rpA loop
        if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag < dist then
          iii := ng - rp - ii + 1;
          dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag;
        end if;
      end for;
      h := gammaReordered[ng - rp - 2*i + 1];
      gammaReordered[ng - rp - 2*i + 1] := gammaReordered[iii];
      gammaReordered[iii] := h;

//###############################################################

      Modelica.Utilities.Streams.print("2x2, compl system, real ass Ann = "+Matrices.printMatrix(A_rsf[n-1:n,n-1:n])+"\n ap = ");

      Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("as", gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)]);

//###############################################################

      Ks2 := Internal.assignOneOrTwoPoles(
        A_rsf[n - 1:n, n - 1:n],
        matrix(B_rsf[n - 1:n, :]),
        gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)],
        tolerance);

      K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
      A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
      select := fill(false, n - counter + 1);
      select[n - counter:n - counter + 1] := {true,true};

      (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
        "E",
        "V",
        select,
        A_rsf[counter:n, counter:n],
        identity(n - counter + 1)); //The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
      B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
      ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
      counter := counter + 2;
      counter2 := counter2 + 2;
    end for;
  end if;

  if counter2>rpg and counter2<rpA then//System has more real eigenvalues than real assigned poles
    for i in 1:div(rpA - rpg, 2) loop// 2x2 blocks; 2 real system poles and a pair of complex assigned poles; take the next two real
                                      // eigenvalues in the diagonal of the Schur form and search the complex pole pair of the assigned poles
                                      // which is nearest to the two real poles
      dist := Modelica.Constants.inf;
      for ii in 2*(i - 1)+1:2:rpA - rpg loop
//        if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[n - 1, n - 1] -
//          gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].re) < dist then
        if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[n - 1, n - 1] -
          gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].im) < dist then
          iii := ng - rp - ii + 1;
          dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re)
             + abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[
            n - 1, n - 1] - gammaReordered[ng - rp - ii + 1].re) +
            abs(gammaReordered[ng - rp - ii + 1].im);
        end if;
      end for;
      h := gammaReordered[ng - rp - 2*(i - 1)];
      gammaReordered[ng - rp - 2*(i - 1)] := gammaReordered[iii];
      gammaReordered[iii] := h;
      h := gammaReordered[ng - rp - 2*i + 1];
      gammaReordered[ng - rp - 2*i + 1] := gammaReordered[iii - 1];
      gammaReordered[iii - 1] := h;

//###############################################################
Modelica.Utilities.Streams.print("2x2, 2 real system, complex ass Ann = "+Matrices.printMatrix(A_rsf[n-1:n,n-1:n])+"\n ap = "+String(h));
//###############################################################

      Ks2 := Internal.assignOneOrTwoPoles(
        A_rsf[n - 1:n, n - 1:n],
        matrix(B_rsf[n - 1:n, :]),
        gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)],
        tolerance);

      K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
      A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
      select := fill(false, n - counter + 1);
      select[n - counter:n - counter + 1] := {true,true};

      (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
        "E",
        "V",
        select,
        A_rsf[counter:n, counter:n],
        identity(n - counter + 1)); //The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
      B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
      ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
      counter := counter + 2;
      counter2 := counter2 + 2;
      Modelica.Utilities.Streams.print("counter2Case3 = " + String(counter2));
    end for;
  end if;
//      else

  for i in 1:ncc loop // 2x2 blocks; 2 complex system poles and two complex assigned poles; take the next complex
                      // system pole pair (next Schur bump) in the diagonal of the Schur form and search the complex
                      //  assigned pole pair which is nearest
    dist := Modelica.Constants.inf;
    evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);//positive imaginary part of the complex system pole pair
    for ii in 2*(i - 1) + 1:2:2*ncc loop
      if abs(A_rsf[n, n] - gammaReordered[2*ncc - ii + 1].re) + abs(evImag -
          abs(gammaReordered[2*ncc - ii + 1].im)) < dist then
        iii := 2*ncc - ii + 1;
        dist := abs(A_rsf[n, n] - gammaReordered[2*ncc - ii + 1].re) + abs(
          evImag - abs(gammaReordered[2*ncc - ii + 1].im));
      end if;
    end for;
    h := gammaReordered[2*ncc - 2*(i - 1)];
    gammaReordered[2*ncc - 2*(i - 1)] := gammaReordered[iii];
    gammaReordered[iii] := h;
    h := gammaReordered[2*ncc - 2*i + 1];
    gammaReordered[2*ncc - 2*i + 1] := gammaReordered[iii - 1];
    gammaReordered[iii - 1] := h;

    //###############################################################
Modelica.Utilities.Streams.print("2x2, 2 compl system, complex ass Ann = "+Matrices.printMatrix(A_rsf[n-1:n,n-1:n])+"\n ap = "+String(h));
//###############################################################

    Ks2 := Internal.assignOneOrTwoPoles(
      A_rsf[n - 1:n, n - 1:n],
      matrix(B_rsf[n - 1:n, :]),
      gammaReordered[2*ncc - 2*i + 1:2*ncc - 2*(i - 1)],
      tolerance);
    K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
    A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
    select := fill(false, n - counter + 1);
    select[n - counter:n - counter + 1] := {true,true};

    (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
      "E",
      "V",
      select,
      A_rsf[counter:n, counter:n],
      identity(n - counter + 1));   //The Schur vector matrix is identity, since A_rsf already has Schur form

    A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
    B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
    ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
    counter := counter + 2;
    counter2 := counter2 + 2;
  end for;

  S := ss.A - ss.B*K;
  po :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(S);

  if calculateEigenvectors then
//     X := fill(Complex(0), n, n);
//     for i in 1:n loop
//       SS := Complex(1)*S;
//       for ii in 1:n loop
//         SS[ii, ii] := SS[ii, ii] - po[i];
//       end for;
//       Xj := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
//                                  SS);
//       for ii in 1:n loop
//         X[ii, i] := Xj[ii, 1];
//       end for;
//     end for;
//      Modelica_LinearSystems2.Math.Complex.Matrices.print(X,6,"X1");
    X :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenVectors(S);
//      Modelica_LinearSystems2.Math.Complex.Matrices.print(X,6,"X2");

  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  (K, S, po, nfp, nap, nup) </td><td align=center> =  </td>  <td> StateSpace.Design.<b>assignPolesMI</b>(ss, gamma, np, tol, calculateEigenvectors)  </td> </tr>
</table>

<h4>Description</h4>
<p>
The purpose of this function is to determine the state feedback matrix <b>K</b> for a
given time invariant multi input state system (<b>A</b>,<b>B</b>) such that the
closed-loop state matrix <b>A</b>-<b>B</b>*<b>K</b> has specified eigenvalues. The
feedback matrix <b>K</b> is calculated by factorization following [1]. The algorithm
modifies the eigenvalues sequentially and also allows partial eigenvalue assignment.<br>
</p>
<p>
At the beginning of the algorithm, the feedback matrix <b>K</b> is set to zero (<b>K</b> = <b>0</b>) and the matrix <b>A</b> is
reduced to an ordered real Schur form by separating its spectrum in two parts
</p>
<blockquote><pre>
             | <b>F</b>1  <b>F</b>3|
<b>F</b> = <b>Q</b>*<b>A</b>*<b>Q</b>' = |       |
             | <b>0</b>   <b>F</b>2|
</pre></blockquote>
<p>
in such a way, that <b>F</b>1 contains the eigenvalues that will be
retained and <b>F</b>3 contains the eigenvalues going to be modified. On the suggestion
of [1] the eigenvalues <i>evr</i> to be retained are chosen as
</p>
<blockquote><pre>
evr = {s in C: Re(s) &lt; -alpha, alpha &gt; =0}
</pre> </blockquote>
<p>
but other specification are conceivable of course.<br>
</p>

<p>
Let
</p>
<blockquote><pre>
<b>G</b> = [<b>G</b>1;<b>G</b>2] = <b>Q</b>*<b>B</b>
</pre> </blockquote>
<p>
with an appropriate partition according to <b>F</b>2. (<b>F</b>2, <b>G</b>2) has to be
controllable.
</p>

<p>
If the feedback matrix <b>K</b> is taken in a form
</p>
<blockquote><pre>
<b>K</b> = [0, <b>K</b>2]
</pre></blockquote>
<p>
the special structure of <b>F</b> and <b>K</b> results in a closed loop state
matrix
</p>
<blockquote><pre>
          |<b>F</b>1 <b>F</b>3 - <b>G</b>1*<b>K</b>2|
<b>F</b> - <b>G</b>*<b>K</b> = |             |
          |0  <b>F</b>2 - <b>G</b>2*<b>K</b>2|
</pre></blockquote>
<p>
with only the eigenvalues of <b>F</b>2 are modified. This approach to modify
separated eigenvalues is used to sequentially shift one real eigenvalue ore two
complex conjugated eigenvalues stepwise until all assigned eigenvalues are placed.
Therefore, at each step i always the (two) lower right eigenvalue(s) are modified by an
appropriate feedback matrix <b>K</b>i. The matrix <b>F</b> - <b>G</b>*<b>K</b>i remains in real Schur form. The
assigned eigenvalue(s) is (are) then moved to another diagonal position of the real Schur
form using reordering techniques <b>F</b> &lt; -- <b>Q</b>i*<b>F</b>*<b>Q</b>i'  and a new block is transferred to the
lower right diagonal position. The transformations are accumulated in <b>Q</b>i and are also
applicated to the matrices
</p>
<blockquote><pre>
<b>G</b> &lt; - <b>Q</b>i*<b>G</b> <b>Q</b> &lt; - <b>Q</b>i*<b>Q</b>
</pre></blockquote>
<p>
The eigenvalue(s) to be assigned at  each step is (are) chosen such that the norm of each <b>K</b>i is minimized [1].
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Real Q[3,3];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.StateSpace.Analysis.observabilityMatrix(ss);
// Q = [0, 1, 0; 0, 1, 1; 1, 1, 2]
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga A. (1981):</dt>
<dd> <b>A Schur method for pole assignment</b>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>

</html>"));
end assignPolesMI;
