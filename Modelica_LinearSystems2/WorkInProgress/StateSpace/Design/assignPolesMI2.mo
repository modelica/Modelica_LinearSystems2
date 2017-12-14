within Modelica_LinearSystems2.WorkInProgress.StateSpace.Design;
encapsulated function assignPolesMI2
  "Pole assignment design algorithm for multi input systems"

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Matrices;
  import Complex;

  input StateSpace ss "state space system";

  input Complex gamma[:]=fill(Complex(0), size(ss.A,1)) "Designed Poles";
  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
  input Boolean giveEigenvalues=false
    "Boolean parameter to display the order of the eigenvalues";
  input Real alpha=-1e10
    "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
  input Real tolerance=Modelica.Math.Matrices.norm(ss.A, 1)*1e-12
    "The tolerance to be used in determining the controllability of (A,B)";
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
  Integer counter;
  Integer n=size(ss.A, 1);
  Integer ipf;
  Integer nccA "number of conjugated complex pole pairs of unmodified system";
  Integer nccg "number of conjugated complex pole pairs of gamma";
  Integer ncc;

  Complex h;
  Integer h2[2];
  Real alphaReal[size(ss.A, 1)] "Real part of eigenvalue=alphaReal+i*alphaImag";
  Real alphaImag[size(ss.A, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

  Complex pi[:]=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);

  Boolean complex_assignedPoles=false
    "True, if there is at least one conjugated comples pole pair in the set of the assigned poles";
  Boolean complex_originalPoles=false
    "True, if there is at least one conjugated comples pole pair in the set of unmodified system poles";
  Boolean consistency;

  Complex SS[:,:];
  Complex Xj[:,:];

  Integer markA[n]=fill(1,n);
  Integer markg[n]=fill(1,n);

  Complex ev[:];

algorithm
  if giveEigenvalues then
    (,,alphaReal,alphaImag,info) :=
      Modelica_LinearSystems2.Math.Matrices.LAPACK.dgees(ss.A);
    assert(info == 0, "The output info of LAPACK.dgees should be zero, else if\n
     info < 0:  if info = -i, the i-th argument of dgees had an illegal value\n
     info > 0:  if INFO = i, and i is
               <= N: the QR algorithm failed to compute all the
                     eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
                     contain those eigenvalues which have converged; if
                     JOBVS = 'V', VS contains the matrix which reduces A
                     to its partially converged Schur form.\n");
    for i in 1:n loop
      po[i].re := alphaReal[i];
      po[i].im := alphaImag[i];
    end for;
    Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("The eigenvalues of the open loop system are sorted to\n eigenvalues", pi);
  else
    assert(size(gamma, 1) <= size(ss.A, 1), "At most n (order of ss) eigenvalues can be assigned");

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

  // determine number of poles not to be assigned according to alpha
     nfp := 0;
     nccA := 0;
     nccg := 0;
     for i in 1:n loop
       if alphaReal[i] < alpha then
         nfp := nfp + 1;
         markA[i] := 0;
         markg[i] := 0;
       end if;
       if abs(alphaImag[i]) > 0 then
         markA[i] := 2;
         nccA := nccA + 1;
       end if;
       if abs(gammaReordered[i].im) > 0 then
         markg[i] := 2;
         nccg := nccg + 1;
       end if;
     end for;
     nap := n - nfp;
     nccA := div(nccA, 2);
     nccg := div(nccg, 2);
     ncc := max(nccA, nccg);

   // reorder real Schur form according to alpha
   (A_rsf,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSFc(
       A_rsf,
       identity(size(A_rsf, 1)),
       alphaReal,
       alphaImag,
       alpha);
   ZT := transpose(Z)*ZT;
   B_rsf := ZT*ss.B;

    Modelica_LinearSystems2.Math.Vectors.printVector(alphaReal,6,"alphaReal");
    Modelica_LinearSystems2.Math.Vectors.printVector(alphaImag,6,"alphaImag");
    Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("gammaReordered1", gammaReordered);

  // Reorder gammaReordered according to alpha
    ii := 1;
    for i in 1:n loop
      if markg[i]==0 then
        h := gammaReordered[ii];
        gammaReordered[ii] := gammaReordered[i];
        gammaReordered[i] := h;
        ii := ii+1;
      end if;
    end for;

   // check consistency of poles assignment, i.e. complex pole pairs in gammaReordered and alphaReal and alphaImag may not be separated
    consistency:=true;
    i:=1;
    while i<n loop
      if markA[i]==2 then
        consistency :=markA[i + 1] == 2 and ((markg[i] == 2 and markg[i + 1] == 2)
          or (markg[i] == 1 and markg[i + 1] == 1));
        i := i+2;
      elseif markg[i]==2 then
        consistency :=markg[i + 1] == 2 and (markA[i] == 1 and markA[i + 1] == 1);
        i := i+2;
      else
        i := i+1;
      end if;
    end while;
   assert(consistency,"System poles and assigned poles have to be assigned consistently, i.e. complex pole pairs may not be separated");

  // main algorithm
  K := zeros(size(ss.B, 2), size(ss.A, 1));
  counter := nfp + 1;
   while counter<=n loop
    if markA[n + nfp + 1 - counter] == 2 or markg[n + nfp + 1 - counter] == 2 then

//#############################
      Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("g2", gammaReordered[n + nfp - counter:n + nfp + 1 - counter]);
      Modelica_LinearSystems2.Math.Vectors.printVector(alphaReal[n + nfp - counter:n + nfp + 1 - counter],6,"ar2");
      ev:=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(A_rsf[n - 1:n, n - 1:n]);
      Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("ev2", ev);
//#############################

      Ks2 := StateSpace.Internal.assignOneOrTwoPoles(
        A_rsf[n - 1:n, n - 1:n],
        matrix(B_rsf[n - 1:n, :]),
        gammaReordered[n + nfp - counter:n + nfp +
        1 - counter],
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
        identity(n - counter + 1));  //The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
      B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
      ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
//       h2 := markA[counter:counter + 1];
//       markA[counter:counter + 1] := markA[n - 1:n];
//       markA[n - 1:n] := h2;
//       h2 := markg[counter:counter + 1];
//       markg[counter:counter + 1] := markg[n - 1:n];
//       markg[n - 1:n] := h2;
//       h := gammaReordered[n - 1];
//       gammaReordered[n - 1] := gammaReordered[counter];
//       gammaReordered[counter] := gammaReordered[n - 1];
//       h := gammaReordered[n];
//       gammaReordered[n] := gammaReordered[counter + 1];
//       gammaReordered[counter + 1] := gammaReordered[n];
      counter := counter + 2;
    else
//#############################
print("g1 = "+String(gammaReordered[n + nfp + 1 - counter]));
print("ar1 = "+String(alphaReal[n + nfp + 1 - counter]));
print("ev1 = "+String(A_rsf[n,n])+"\n");
//#############################

      Ks1 := StateSpace.Internal.assignOneOrTwoPoles(
        matrix(A_rsf[n, n]),
        transpose(matrix(B_rsf[n, :])),
        {gammaReordered[n + nfp + 1 - counter]},
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
        identity(n - counter + 1)); //The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q1;
      B_rsf[counter:n, :] := transpose(Q1)*B_rsf[counter:n, :];
      ZT[counter:n, :] := transpose(Q1)*ZT[counter:n, :];
//       h2[1] := markA[counter];
//       markA[counter] := markA[n];
//       markA[n] := h2[1];
//       h2[1] := markg[counter];
//       markg[counter] := markg[n];
//       markg[n] := h2[1];
//       h := gammaReordered[n];
//       gammaReordered[n] := gammaReordered[counter];
//       gammaReordered[counter] := gammaReordered[n];
      counter := counter + 1;
    end if;

    end while;

  S := ss.A - ss.B*K;
  po :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(S);
   X :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenVectors(S);

//    X := fill(Complex(0),n,n);
//    for i in 1:n loop
//      SS:=Complex(1)*S;
//      for ii in 1:n loop
//        SS[ii,ii] := SS[ii,ii]-po[i];
//      end for;
//      Xj := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
//                                 SS);
//      // Matrices.printMatrix(Complex.real(Xj),6,"ReXj");
//      // Matrices.printMatrix(Complex.imag(Xj),6,"ImXj");
//      for ii in 1:n loop
//        X[ii,i] := Xj[ii,1];
//      end for;
//    end for;

end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  (K, S, po, nfp, nap, nup) </td><td align=center> =  </td>  <td> StateSpace.Design.<b>assignPolesMI</b>(ss, gamma, np, alpha, tol)  </td> </tr>
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
At the beginning of the algorithm, the feedback matrix <b>K</b> is set to zero
(<b>K</b> = <b>0</b>) and the matrix <b>A</b> is reduced to an ordered real Schur
form by separating its spectrum in two parts
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
Therefore, at each step i always the (two) lower right eingenvalue(s) are modified by an
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
end assignPolesMI2;
