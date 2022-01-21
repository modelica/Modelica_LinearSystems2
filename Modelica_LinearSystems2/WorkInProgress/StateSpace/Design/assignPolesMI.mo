within Modelica_LinearSystems2.WorkInProgress.StateSpace.Design;
encapsulated function assignPolesMI
  "Pole assignment design algorithm for multi input systems"

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import MatricesMSL = Modelica.Math.Matrices;
  import Complex;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.ComplexMathAdds;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;

  input StateSpace ss "state space system";

  input Complex gamma[:]=fill(Complex(0), 0) "Designed Poles";
//  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
  input Real alpha=-1e10
    "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
  input Real tolerance=MatricesMSL.norm(ss.A, 1)*1e-12
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
  (A_rsf,Z,alphaReal,alphaImag) := MatricesMSL.realSchur(ss.A);
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
  gammaReordered := Modelica.ComplexMath.Vectors.reverse(gammaReordered);
  nccg := div(size(gammaReordered, 1) - rpg, 2);
  ncc := min(nccA, nccg);
  rp := min(rpA, rpg);
  if nccA > 0 then
    (A_rsf[nfp + 1:n, nfp + 1:n],Q2) := MatricesMSL.LAPACK.dtrsen(
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

      (A_rsf[counter:n, counter:n],Q1) := MatricesMSL.LAPACK.dtrsen(
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

  ComplexMathAdds.Vectors.print("as",gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)]);

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

      (A_rsf[counter:n, counter:n],Q2) := MatricesMSL.LAPACK.dtrsen(
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

      (A_rsf[counter:n, counter:n],Q2) := MatricesMSL.LAPACK.dtrsen(
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

    (A_rsf[counter:n, counter:n],Q2) := MatricesMSL.LAPACK.dtrsen(
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
  po := ComplexMathAdds.eigenValues(S);

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
//      ComplexMathAdds.Matrices.print(X,6,"X1");
    X := ComplexMathAdds.eigenVectors(S);
//      ComplexMathAdds.Matrices.print(X,6,"X2");

  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  (K, S, po, nfp, nap, nup) </td><td align=center> =  </td>  <td> StateSpace.Design.<strong>assignPolesMI</strong>(ss, gamma, np, tol, calculateEigenvectors)  </td> </tr>
</table>

<h4>Description</h4>
<p>
The purpose of this function is to determine the state feedback matrix <strong>K</strong> for a
given time invariant multi input state system (<strong>A</strong>,<strong>B</strong>) such that the
closed-loop state matrix <strong>A</strong>-<strong>B</strong>*<strong>K</strong> has specified eigenvalues. The
feedback matrix <strong>K</strong> is calculated by factorization following [1]. The algorithm
modifies the eigenvalues sequentially and also allows partial eigenvalue assignment.<br>
</p>
<p>
At the beginning of the algorithm, the feedback matrix <strong>K</strong> is set to zero (<strong>K</strong> = <strong>0</strong>) and the matrix <strong>A</strong> is
reduced to an ordered real Schur form by separating its spectrum in two parts
</p>
<blockquote><pre>
             | <strong>F</strong>1  <strong>F</strong>3|
<strong>F</strong> = <strong>Q</strong>*<strong>A</strong>*<strong>Q</strong>' = |       |
             | <strong>0</strong>   <strong>F</strong>2|
</pre></blockquote>
<p>
in such a way, that <strong>F</strong>1 contains the eigenvalues that will be
retained and <strong>F</strong>3 contains the eigenvalues going to be modified. On the suggestion
of [1] the eigenvalues <em>evr</em> to be retained are chosen as
</p>
<blockquote><pre>
evr = {s in C: Re(s) &lt; -alpha, alpha &gt; =0}
</pre></blockquote>
<p>
but other specification are conceivable of course.<br>
</p>

<p>
Let
</p>
<blockquote><pre>
<strong>G</strong> = [<strong>G</strong>1;<strong>G</strong>2] = <strong>Q</strong>*<strong>B</strong>
</pre></blockquote>
<p>
with an appropriate partition according to <strong>F</strong>2. (<strong>F</strong>2, <strong>G</strong>2) has to be
controllable.
</p>

<p>
If the feedback matrix <strong>K</strong> is taken in a form
</p>
<blockquote><pre>
<strong>K</strong> = [0, <strong>K</strong>2]
</pre></blockquote>
<p>
the special structure of <strong>F</strong> and <strong>K</strong> results in a closed loop state
matrix
</p>
<blockquote><pre>
          |<strong>F</strong>1 <strong>F</strong>3 - <strong>G</strong>1*<strong>K</strong>2|
<strong>F</strong> - <strong>G</strong>*<strong>K</strong> = |             |
          |0  <strong>F</strong>2 - <strong>G</strong>2*<strong>K</strong>2|
</pre></blockquote>
<p>
with only the eigenvalues of <strong>F</strong>2 are modified. This approach to modify
separated eigenvalues is used to sequentially shift one real eigenvalue ore two
complex conjugated eigenvalues stepwise until all assigned eigenvalues are placed.
Therefore, at each step i always the (two) lower right eigenvalue(s) are modified by an
appropriate feedback matrix <strong>K</strong>i. The matrix <strong>F</strong> - <strong>G</strong>*<strong>K</strong>i remains in real Schur form. The
assigned eigenvalue(s) is (are) then moved to another diagonal position of the real Schur
form using reordering techniques <strong>F</strong> &lt; -- <strong>Q</strong>i*<strong>F</strong>*<strong>Q</strong>i'  and a new block is transferred to the
lower right diagonal position. The transformations are accumulated in <strong>Q</strong>i and are also
applicated to the matrices
</p>
<blockquote><pre>
<strong>G</strong> &lt; - <strong>Q</strong>i*<strong>G</strong> <strong>Q</strong> &lt; - <strong>Q</strong>i*<strong>Q</strong>
</pre></blockquote>
<p>
The eigenvalue(s) to be assigned at  each step is (are) chosen such that the norm of each <strong>K</strong>i is minimized [1].
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Real Q[3,3];

<strong>algorithm</strong>
  Q := Modelica_LinearSystems2.StateSpace.Analysis.observabilityMatrix(ss);
// Q = [0, 1, 0; 0, 1, 1; 1, 1, 2]
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga A. (1981):</dt>
<dd> <strong>A Schur method for pole assignment</strong>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>

</html>"));
end assignPolesMI;
