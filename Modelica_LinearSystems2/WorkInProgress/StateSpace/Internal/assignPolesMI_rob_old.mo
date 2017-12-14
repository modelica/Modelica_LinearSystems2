within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function assignPolesMI_rob_old
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Complex;
  import matMul = Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.matMatMul;
  import Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.matVecMul;
  import Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.C_transpose;
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;
  import Modelica.Utilities.Streams.print;

  input Real A[:,size(A, 1)] "system matrix";
  input Real B[size(A, 1),:] "control input matrix";
  input Complex gamma[size(A, 1)];
  input Boolean IniX=false "Initial values of X are provided";
  input Complex Xini[size(A, 1),size(A, 1)]=fill(Complex(0),size(A, 1),size(A, 1))
    "Initial values of the eigenvectors X";

  output Real K[size(B, 2),size(A, 1)] "feedback matrix";
  output Complex evX[:,:] "eigen vectors of the closed loop system";
protected
  Complex AC[size(A, 1),size(A, 2)]=Complex(1)*A;
  Complex Lambda[size(A, 1),size(A, 1)];
  Real Ur[:,:];
  Real Vr[:,:];
  Complex U0[:,:];
  Complex U1[:,:];
  Complex U1T[:,:];
  Complex X[size(A, 1),size(A, 2)]=Xini;
  Complex Xj[size(A, 1),size(A, 1)-1];
  Complex XC[size(A, 1),size(A, 2)];
  Complex XC2[size(A, 1),size(A, 2)];
  Complex Z[:,:];
  Complex M[size(A, 1),size(A, 1)];
  Complex KC[:,:];
  Complex QX[:,:];
  Complex C[:,:];
  Complex Cc[:,:];
  Complex Sr[:,:];
  Complex ST[:,:];

  Complex gammaSorted[size(gamma, 1)];
  Complex gammaSorted2[size(gamma, 1)];
  Real sigmaB[:];
  Complex ev[size(A, 1)];

  Complex y[:];

  Real condX2;
  Real norm_y;
  Integer rankB;
  Integer nx=size(A, 1);
  Integer numberOfRealEigenvalues;
  Integer numberOfComplexPairs;
  Integer i;
  Integer l1;
  Integer l2;
  Integer k;
  Integer idx;
  Real eps=Modelica.Constants.eps;
  Integer maxSteps=3;

  Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal.assignPolesMI_rob_old.subSpace
    subS[                                                                         size(gamma,1)];

  Complex a;
  Complex MM[:,:];
  Integer cnt;
  StateSpace ss=StateSpace(A=A, B=B, C=zeros(1,nx), D=zeros(1,size(B,2)));

algorithm
  //check controllability
  assert(StateSpace.Analysis.isControllable(ss),"Poles cannot be placed since system is not controllable");

  // sort eigenvalues to [real ev, complex ev(im>0), conj(complex ev(im>0))]
  (gammaSorted2,numberOfRealEigenvalues) :=
    Modelica_LinearSystems2.Internal.reorderZeros(gamma);
  gammaSorted := gammaSorted2;
  numberOfComplexPairs := integer((nx - numberOfRealEigenvalues)/2);
  for i in numberOfRealEigenvalues + 1:numberOfRealEigenvalues + numberOfComplexPairs loop
    gammaSorted[i] :=if Im(gammaSorted2[2*i - numberOfRealEigenvalues - 1]) > 0 then gammaSorted2[2*i - numberOfRealEigenvalues - 1] else Modelica.ComplexMath.conj(gammaSorted2[2*i - numberOfRealEigenvalues - 1]);
    gammaSorted[i + numberOfComplexPairs] :=Modelica.ComplexMath.conj(gammaSorted[i]);
  end for;

  for i in 1:nx loop
    Lambda[i, i] := gammaSorted[i];
  end for;

//   (sigmaB,Q,R) := Matrices.C_singularValues(BC);
//   rankB := 0;
//   i := size(sigmaB, 1);
//   while i > 0 loop
//     if sigmaB[i] > 1e-10 then
//       rankB := i;
//       i := 0;
//     end if;
//     i := i - 1;
//   end while;

  (sigmaB,Ur,Vr) := Modelica.Math.Matrices.singularValues(B);
  rankB := 0;
  i := size(sigmaB, 1);
  while i > 0 loop
    if sigmaB[i] > 1e-10 then
      rankB := i;
      i := 0;
    end if;
    i := i - 1;
  end while;

//    R := C_transpose(R);
//    Z := R[1:rankB, 1:size(B, 2)];
//   for l1 in 1:rankB loop
//     for l2 in 1:size(B, 2) loop
//       Z[l1, l2] := Z[l1, l2]/sigmaB[l2];
//     end for;
//   end for;

   Z := fill(Complex(0), size(B, 2), rankB);
   for l1 in 1:rankB loop
     for l2 in 1:size(B, 2) loop
//       Z[l1, l2] := Modelica.ComplexMath.conj(R[l2, l1])/sigmaB[l2];
       Z[l1, l2] := Complex((Vr[l2, l1])/sigmaB[l2]);
     end for;
   end for;

//   U0 := Q[:, 1:rankB];
//   U1 := Q[:, rankB + 1:nx];
//   U1T := C_transpose(U1);

  U0 := Complex(1)*Ur[:, 1:rankB];
  U1 := Complex(1)*Ur[:, rankB + 1:nx];
  U1T := Complex(1)*transpose(Ur[:, rankB + 1:nx]);

  condX2 := eps + 1;

  if numberOfComplexPairs > 0 then
//      (sigmaM,Qu,Ru) := Matrices.C_singularValues([U1T;  matMul(U1T,Complex(1)*A)]);
//       rank_ := 0;
//          i := size(sigmaM,1);
//        while i > 0 loop
//          if sigmaM[i] > 1e-10 then
//            rank_ := i;
//            i := 0;
//          end if;
//          i := i - 1;
//        end while;
//      Ru := C_transpose(Ru);
//      Sr := Ru[:,rank_+1:nx];

     Sr :=  Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                                 [U1T;  matMul(U1T,Complex(1)*A)]);
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(Sr),6,"ReSr");
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(Sr),6,"ImSr");
  else
    Sr := fill(Complex(0),nx,0);
  end if;

  for l1 in 1:nx - numberOfComplexPairs loop

    AC := Complex(-1)*A;

    for l2 in 1:nx loop
      AC[l2, l2] := AC[l2, l2] + gammaSorted[l1];
    end for;

    C := matMul(U1T, AC);

    Cc := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                               [C; C_transpose(Sr)]);

//     if numberOfComplexPairs > 0 and l1 > numberOfRealEigenvalues then
// //       (sigmaC,QC,RC) := Matrices.C_singularValues([C; C_transpose(Sr)]);
// //       rankC := 0;
// //       i := size(sigmaC, 1);
// //       while i > 0 loop
// //         if sigmaC[i] > 1e-10 then
// //           rankC := i;
// //           i := 0;
// //         end if;
// //         i := i - 1;
// //       end while;
// //       RC := C_transpose(RC);
// //       subS[l1].S := [RC[:, rankC + 1:nx],Sr];
//
//       RC := Matrices.C_nullspace([C; C_transpose(Sr)]);
//       subS[l1].S := [RC,Sr];
//     else
// //       (sigmaC,QC,RC) := Modelica_LinearSystems2.Math.Matrices.C_singularValues(C);
// //       rankC := 0;
// //       i := size(sigmaC, 1);
// //       while i > 0 loop
// //         if sigmaC[i] > 1e-10 then
// //           rankC := i;
// //           i := 0;
// //         end if;
// //         i := i - 1;
// //       end while;
// //       RC := C_transpose(RC);
// //       subS[l1].S := RC[:, rankC + 1:nx];
//
//       subS[l1].S := Matrices.C_nullspace(C);
//     end if;

    subS[l1].S := if l1 > numberOfRealEigenvalues then [Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                                                                             [C; C_transpose(Sr)]),Sr] else Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                                                                                            C);

//    subS[l1].S := Matrices.C_nullspace(C);
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(subS[l1].S),6,"ResubS[l1].S");
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(subS[l1].S),6,"ImsubS[l1].S");

// MM := matMul(C,Modelica.ComplexMath.conj(Matrices.C_nullspace(Cc)));
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(MM),6,"ReMM");
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(MM),6,"ImMM");

// MM := matMul(C,Sr);
// //MM := matMul(C_transpose(Sr),Cc);
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Re(MM),6,"ReMM");
// Modelica_LinearSystems2.Math.Matrices.printMatrix(Im(MM),6,"ImMM");
end for;

if not IniX then
// initialization of X according to Byers
  for l1 in 1:nx - numberOfComplexPairs loop
    y := fill(Complex(0), nx);
    for l2 in 1:size(subS[l1].S, 2) loop
      y := X[:, l1] + X[:, l1] + subS[l1].S[:, l2];
    end for;
    y :=y/Modelica.ComplexMath.Vectors.norm(y);
    for l2 in 1:nx loop
      X[l2, l1] := y[l2];
    end for;
  end for;
  for l1 in 1:numberOfComplexPairs loop
    for l2 in 1:nx loop
      X[l2, numberOfRealEigenvalues + numberOfComplexPairs + l1] :=Modelica.ComplexMath.conj(X[l2, numberOfRealEigenvalues + l1]);
    end for;
  end for;

 // initialization of X according to place.m
   for l1 in 1:nx - numberOfComplexPairs loop
     for l2 in 1:nx loop
       X[l2, l1] := subS[l1].S[l2, 1];
     end for;
   end for;
   for l1 in 1:numberOfComplexPairs loop
     for l2 in 1:nx loop
       X[l2, numberOfRealEigenvalues + numberOfComplexPairs + l1] :=Modelica.ComplexMath.conj(X[l2, numberOfRealEigenvalues + l1]);
     end for;
   end for;
end if;
// eigenvector modification
  k := 0;
  while (k < maxSteps) loop
    k := k + 1;

    for l1 in 1:nx - numberOfComplexPairs loop
      if l1 == 1 then
        for l2 in 1:nx loop
          for l3 in 1:nx - 1 loop
            Xj[l2, l3] := X[l2, l3 + 1];
          end for;
        end for;
      else
        for l2 in 1:nx loop
          for l3 in 1:l1 - 1 loop
            Xj[l2, l3] := X[l2, l3];
          end for;
          for l3 in l1:nx - 1 loop
            Xj[l2, l3] := X[l2, l3 + 1];
          end for;
        end for;
      end if;

      QX := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_QR2(Xj);
      ST := C_transpose(subS[l1].S);
      y := matVecMul(ST, QX[:, nx]);
      norm_y :=Modelica.ComplexMath.Vectors.norm(y);
      y := matVecMul(subS[l1].S, y)/norm_y;

//        if l1 > numberOfRealEigenvalues and Complex.'abs'(Complex.Vectors.multiply(y,Modelica.ComplexMath.conj(y)))>0.9 then
//          idx := 1 + rem(k, size(subS[l1].S, 2) - size(Sr, 2));
//          print(" k = "+String(k)+", l1 = "+String(l1)+", idx = " + String(idx));
//          y := (y + subS[l1].S[:, idx])/sqrt(2);
//        end if;

      for l2 in 1:nx loop
        X[l2, l1] := y[l2];
      end for;

      if l1 > numberOfRealEigenvalues then
        for l2 in 1:nx loop
          X[l2, l1 + numberOfComplexPairs] :=Modelica.ComplexMath.conj(y[l2]);
        end for;
      end if;

        end for;
condX2 := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.conditionNumber(
                                           X);

//print("\ncondX2 = "+String(condX2));
  end while;

  XC := C_transpose(X);
  XC2 := C_transpose(matMul(X, Lambda));
  M := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_solve2(
                                                      XC, XC2);
  M := C_transpose(M);

  for l2 in 1:nx loop
    for l3 in 1:nx loop
      M[l2, l3] := M[l2, l3] - A[l2, l3];
    end for;
  end for;

  KC := matMul(Z, matMul(C_transpose(U0), M));
  K := -Re(KC);
  evX := X;

  ev := Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(A - B*K);
//    Complex.Vectors.print("gammaSorted", gammaSorted);
//    Complex.Vectors.print("ev", ev);

public
  encapsulated record subSpace
    import Modelica;
    import Modelica_LinearSystems2;
    import Complex;
    extends Modelica.Icons.Record;
    Complex S[:,:];
  end subSpace;

end assignPolesMI_rob_old;
