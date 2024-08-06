within Modelica_LinearSystems2.WorkInProgress.StateSpace.Design;
encapsulated function assignPolesMI2
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

  input StateSpace ss "state space system";

  input Complex gamma[:]=fill(Complex(0), size(ss.A,1)) "Designed Poles";
  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
  input Boolean giveEigenvalues=false
    "Boolean parameter to display the order of the eigenvalues";
  input Real alpha=-1e10
    "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
  input Real tolerance=MatricesMSL.norm(ss.A, 1)*1e-12
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

  Complex pi[:]=ComplexMathAdds.eigenValues(ss.A);

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
      MatricesMSL.LAPACK.dgees(ss.A);
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
    ComplexMathAdds.Vectors.print("The eigenvalues of the open loop system are sorted to\n eigenvalues", pi);
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
  (A_rsf,Z,alphaReal,alphaImag) := MatricesMSL.realSchur(ss.A);
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

   Modelica.Math.Vectors.toString(alphaReal, "alphaReal", 6);
   Modelica.Math.Vectors.toString(alphaImag, "alphaImag", 6);
   ComplexMathAdds.Vectors.print("gammaReordered1",gammaReordered);

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
      ComplexMathAdds.Vectors.print("g2",gammaReordered[n + nfp - counter:n + nfp + 1 - counter]);
        Modelica.Math.Vectors.toString(
          alphaReal[n + nfp - counter:n + nfp + 1 - counter], "ar2", 6);
      ev:=ComplexMathAdds.eigenValues( A_rsf[n - 1:n, n - 1:n]);
      ComplexMathAdds.Vectors.print("ev2",ev);
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

      (A_rsf[counter:n, counter:n],Q2) := MatricesMSL.LAPACK.dtrsen(
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

      (A_rsf[counter:n, counter:n],Q1) := MatricesMSL.LAPACK.dtrsen(
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
  po := ComplexMathAdds.eigenValues(S);
   X := ComplexMathAdds.eigenVectors(S);

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
<blockquote><pre>
(K, S, po, nfp, nap, nup) = StateSpace.Design.<strong>assignPolesMI</strong>(ss, gamma, np, alpha, tol)
</pre></blockquote>

<h4>Description</h4>
<p>
The purpose of this function is to determine the state feedback matrix <strong>K</strong> for a
given time invariant multi input state system (<strong>A</strong>,<strong>B</strong>) such that the
closed-loop state matrix <strong>A</strong>-<strong>B</strong>*<strong>K</strong> has specified eigenvalues. The
feedback matrix <strong>K</strong> is calculated by factorization following [1]. The algorithm
modifies the eigenvalues sequentially and also allows partial eigenvalue assignment.<br>
</p>
<p>
At the beginning of the algorithm, the feedback matrix <strong>K</strong> is set to zero
(<strong>K</strong> = <strong>0</strong>) and the matrix <strong>A</strong> is reduced to an ordered real Schur
form by separating its spectrum in two parts
</p>
<blockquote>
  <table style=\"border:0\">
    <tr>
      <td>   </td>
      <td> | </td>
      <td style=\"text-align:center;\"> <strong>F</strong>1 </td>
      <td> &ensp; </td>
      <td style=\"text-align:center;\"> <strong>F</strong>3 </td>
      <td> | </td>
    </tr>
    <tr>
      <td> <strong>F</strong> = <strong>Q</strong>*<strong>A</strong>*<strong>Q</strong>' = </td>
      <td> | </td>
      <td>   </td>
      <td>   </td>
      <td>   </td>
      <td> | </td>
    </tr>
    <tr>
      <td>   </td>
      <td> | </td>
      <td style=\"text-align:center;\"> <strong>0</strong> </td>
      <td>   </td>
      <td style=\"text-align:center;\"> <strong>F</strong>2 </td>
      <td> | </td>
    </tr>
  </table>
</blockquote>
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
<blockquote>
<strong>G</strong> = [<strong>G</strong>1;<strong>G</strong>2] = <strong>Q</strong>*<strong>B</strong>
</blockquote>
<p>
with an appropriate partition according to <strong>F</strong>2. (<strong>F</strong>2, <strong>G</strong>2) has to be
controllable.
</p>
<p>
If the feedback matrix <strong>K</strong> is taken in a form
</p>
<blockquote>
<strong>K</strong> = [<strong>0</strong>, <strong>K</strong>2]
</blockquote>
<p>
the special structure of <strong>F</strong> and <strong>K</strong> results in a closed loop state
matrix
</p>
<blockquote>
  <table style=\"border:0\">
    <tr>
      <td>   </td>
      <td> |  </td>
      <td style=\"text-align:center;\"> <strong>F</strong>1 </td>
      <td> &ensp; </td>
      <td style=\"text-align:right;\"> <strong>F</strong>3 &minus; <strong>G</strong>1*<strong>K</strong>2 | </td>
    </tr>
    <tr>
      <td> <strong>F</strong> &minus; <strong>G</strong>*<strong>K</strong> = </td>
      <td> | </td>
      <td>   </td>
      <td>   </td>
      <td style=\"text-align:right;\"> | </td>
    </tr>
    <tr>
      <td> </td>
      <td> | </td>
      <td style=\"text-align:center;\"> <strong>0</strong> </td>
      <td>   </td>
      <td style=\"text-align:right;\"> <strong>F</strong>2 &minus; <strong>G</strong>2*<strong>K</strong>2 | </td>
    </tr>
  </table>
</blockquote>
<p>
with only the eigenvalues of <strong>F</strong>2 are modified. This approach to modify
separated eigenvalues is used to sequentially shift one real eigenvalue or two
complex conjugated eigenvalues stepwise until all assigned eigenvalues are placed.
Therefore, at each step i always the (two) lower right eingenvalue(s) are modified by an
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
end assignPolesMI2;
