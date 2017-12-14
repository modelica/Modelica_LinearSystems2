within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function assignPolesMI_rob2
  "Modified KNV algorithm. Works like MATLAB's place.m"
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
  input Complex Xini[size(A, 1),size(A, 1)]=fill(Complex(0), size(A, 1), size(A, 1))
    "Initial values of the eigenvectors X";

  output Real K[size(B, 2),size(A, 1)] "feedback matrix";
  output Complex X[size(A, 1),size(A, 2)]=Xini
    "eigen vectors of the closed loop system";
protected
  Real U0[:,:];
  Real Z[:,:];
  Complex S[:,:];
  Integer maxSteps=5;
  Integer numberOfComplexPairs;
  Integer rankB;

  Complex AC[size(A, 1),size(A, 2)]=Complex(1)*A;
  Complex Lambda[size(A, 1),size(A, 1)];
  Real Ur[:,:];
  Real Vr[:,:];

  Complex U1T[:,:];
  Complex Xjj[size(A, 1),size(A, 1)];
  Complex XC[size(A, 1),size(A, 2)];
  Complex XC2[size(A, 1),size(A, 2)];
  Complex MM[size(A, 1),size(A, 1)];
  Real M[size(A, 1),size(A, 1)];
  Complex KC[:,:];
  Complex QX[:,:];
  Complex C[:,:];
  Complex Cc[:,:];
  Complex Sr[:,:];
  Real Srr[:,:];
  Real U1Tr[:,:];
  Complex ST[:,:];

  Complex S2[size(A, 1),size(B, 2)];

  Complex gammaSorted[size(gamma, 1)];
  Complex gammaSorted2[size(gamma, 1)];
  Real sigmaB[:];
//  Complex ev[size(A, 1)];

  Complex y[:];
  Real condX2;
  Real norm_y;

  Integer nx=size(A, 1);
  Integer numberOfRealEigenvalues;

  Integer i;
  Integer l1;
  Integer l2;
  Integer k;
  Integer idx;
  Real eps=Modelica.Constants.eps;

//  Modelica_LinearSystems2.StateSpace.Internal.assignPolesMI_rob.subSpace subS[size(gamma,1)];

//  Complex a;
//  Integer cnt;
  StateSpace ss=StateSpace(
      A=A,
      B=B,
      C=zeros(1, nx),
      D=zeros(1, size(B, 2)));

   Complex XX[size(A, 1),size(A, 2)];

algorithm
  //check controllability
  assert(StateSpace.Analysis.isControllable(ss), "Poles cannot be placed since system is not controllable");

  // sort eigenvalues to [real ev, complex ev(im>0), conj(complex ev(im>0))]
  (gammaSorted2,numberOfRealEigenvalues) :=
    Modelica_LinearSystems2.Internal.reorderZeros(gamma);
  gammaSorted := gammaSorted2;
  numberOfComplexPairs := integer((nx - numberOfRealEigenvalues)/2);
  for i in numberOfRealEigenvalues + 1:numberOfRealEigenvalues +
      numberOfComplexPairs loop
    gammaSorted[i] :=if Im(gammaSorted2[2*i - numberOfRealEigenvalues - 1]) > 0 then gammaSorted2[2*i - numberOfRealEigenvalues - 1] else Modelica.ComplexMath.conj(gammaSorted2[2*i - numberOfRealEigenvalues - 1]);
    gammaSorted[i + numberOfComplexPairs] :=Modelica.ComplexMath.conj(gammaSorted[i]);
  end for;

  for i in 1:nx loop
    Lambda[i, i] := gammaSorted[i];
  end for;

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

  // ############  Calculation of S (bases of X)  begin
   Z := fill(0, size(B, 2), rankB);
   for l1 in 1:rankB loop
     for l2 in 1:size(B, 2) loop
       Z[l1, l2] := Vr[l2, l1]/sigmaB[l2];
     end for;
   end for;
   U0 := Ur[:, 1:rankB];

 if nx>size(B,2) then

   U1T := Complex(1)*transpose(Ur[:, rankB + 1:nx]);
   U1Tr := transpose(Ur[:, rankB + 1:nx]);
   condX2 := eps + 1;

  S := fill(Complex(0), nx, size(B, 2)*(nx - numberOfComplexPairs));

   if numberOfComplexPairs > 0 and 2*rankB - nx > 0 then
     Srr := Matrices.nullspace([U1Tr; U1Tr*A]);
   else
     Srr := fill(0, nx, 0);
   end if;
   Sr := Complex(1)*Srr;

   //Computation of the nullspaces, i.e. the bases of the eigenvectors
   AC := Complex(-1)*A;
   for l1 in 1:nx - numberOfComplexPairs loop
     for l2 in 1:nx loop
       AC[l2, l2] := Complex(-A[l2, l2]) + gammaSorted[l1];
     end for;

     C := matMul(U1T, AC);

     S2 := if l1 > numberOfRealEigenvalues and 2*rankB - nx > 0 then [Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                                                                                           [C; C_transpose(Sr)]),Sr] else Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
                                                                                            C);
     for l2 in 1:nx loop
       for l3 in 1:rankB loop
         S[l2, rankB*(l1 - 1) + l3] := S2[l2, l3];
       end for;
     end for;
   end for;
 // ############  Calculation of S (bases of X)  end
// (U0, Z, S, rankB) := Modelica_LinearSystems2.StateSpace.Internal.xBase(A, B, gammaSorted, numberOfComplexPairs);

  if not IniX and rankB>1 then
// // initialization of X according to Byers
//   for l1 in 1:nx - numberOfComplexPairs loop
//     y := fill(Complex(0), nx);
//     for l2 in 1:rankB loop
// //      y := X[:, l1] + X[:, l1] + subS[l1].S[:, l2];
//       y := X[:, l1] + X[:, l1] + S[:, rankB*(l1-1)+l2];
//     end for;
//     y := y/Complex.Vectors.norm(y);
//     for l2 in 1:nx loop
//       X[l2, l1] := y[l2];
//     end for;
//   end for;
//   for l1 in 1:numberOfComplexPairs loop
//     for l2 in 1:nx loop
//       X[l2, numberOfRealEigenvalues + numberOfComplexPairs + l1] :=
//         Modelica.ComplexMath.conj(X[l2, numberOfRealEigenvalues + l1]);
//     end for;
//   end for;

 // initialization of X according to place.m
    for l1 in 1:nx - numberOfComplexPairs loop
      for l2 in 1:nx loop
        X[l2, l1] := S[l2, rankB*(l1 - 1) + 1];
      end for;
    end for;
    for l1 in 1:numberOfComplexPairs loop
      for l2 in 1:nx loop
        X[l2, numberOfRealEigenvalues + numberOfComplexPairs + l1] :=Modelica.ComplexMath.conj(X[l2, numberOfRealEigenvalues + l1]);
      end for;
    end for;
  end if;

  if rankB==1 then //X=S
    for l1 in 1:nx - numberOfComplexPairs loop
      for l2 in 1:nx loop
        X[l2, l1] := S[l2, rankB*(l1 - 1) + 1];
      end for;
    end for;
    for l1 in 1:numberOfComplexPairs loop
      for l2 in 1:nx loop
        X[l2, numberOfRealEigenvalues + numberOfComplexPairs + l1] :=Modelica.ComplexMath.conj(X[l2, numberOfRealEigenvalues + l1]);
      end for;
    end for;
  end if;

// eigenvector modification
//  X := Modelica_LinearSystems2.StateSpace.Internal.modifyX( X,S,rankB,numberOfComplexPairs,maxSteps,IniX);

  if rankB>1 then
   k := 0;
   while (k < maxSteps) loop
     k := k + 1;

     for l1 in 1:nx - numberOfComplexPairs loop
       if l1 == 1 then
         for l2 in 1:nx loop
           for l3 in 1:nx - 1 loop
             Xjj[l2, l3] := X[l2, l3 + 1];
           end for;
         end for;
       else
         for l2 in 1:nx loop
           for l3 in 1:l1 - 1 loop
             Xjj[l2, l3] := X[l2, l3];
           end for;
           for l3 in l1:nx - 1 loop
             Xjj[l2, l3] := X[l2, l3 + 1];
           end for;
         end for;
       end if;

       QX := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_QR(
                                                        Xjj);

       ST := C_transpose(S[:, rankB*(l1 - 1) + 1:rankB*l1]);
       y := matVecMul(ST, QX[:, nx]);

       norm_y :=Modelica.ComplexMath.Vectors.norm(y);
       y := matVecMul(S[:, rankB*(l1 - 1) + 1:rankB*l1], y)/norm_y;

          if l1 > numberOfRealEigenvalues and Modelica.ComplexMath.'abs'(Complex.'*'.scalarProduct(y, Modelica.ComplexMath.conj(y))) > 0.9 then
         idx := 1 + rem(k, rankB - size(Sr, 2));
         y := (y + S[:, (l1 - 1)*rankB + idx])/sqrt(2);
       end if;

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
      end while;
    end if;//rankB>1

// Computation of the feedback matrix K
    XC := C_transpose(X);
    XC2 := C_transpose(matMul(X, Lambda));
    MM := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_solve2(
                                                        XC, XC2);
    M := Re(MM);
    M := transpose(M);

    for l2 in 1:nx loop
      for l3 in 1:nx loop
        M[l2, l3] := M[l2, l3] - A[l2, l3];
      end for;
    end for;

  else
    gammaSorted2 :=  Modelica_LinearSystems2.Internal.reorderZeros(gamma);
    for i in 1:nx loop
      Lambda[i, i] := gammaSorted2[i];
    end for;
    X:=fill(Complex(0),nx,nx);

    for i in 1:numberOfRealEigenvalues loop
      M[i,i]:=Re(gammaSorted[i]);
      X[i,i] := Complex(1);
    end for;
    for i in 1:numberOfComplexPairs loop
         M[numberOfRealEigenvalues+2*i-1,numberOfRealEigenvalues+2*i] := -Im(gammaSorted[numberOfRealEigenvalues + 2*i - 1]);
         M[numberOfRealEigenvalues+2*i,numberOfRealEigenvalues+2*i-1] := Im(gammaSorted[numberOfRealEigenvalues + 2*i - 1]);
         M[numberOfRealEigenvalues+2*i-1,numberOfRealEigenvalues+2*i-1] := Re(gammaSorted[numberOfRealEigenvalues + 2*i - 1]);
         M[numberOfRealEigenvalues+2*i,numberOfRealEigenvalues+2*i] := Re(gammaSorted[numberOfRealEigenvalues + 2*i - 1]);
         X[numberOfRealEigenvalues+2*i-1,numberOfRealEigenvalues+2*i-1] := Complex(0.5);
         X[numberOfRealEigenvalues+2*i,numberOfRealEigenvalues+2*i] := Complex(0,0.5);
         X[numberOfRealEigenvalues+2*i-1,numberOfRealEigenvalues+2*i] := Complex(0,-0.5);
         X[numberOfRealEigenvalues+2*i,numberOfRealEigenvalues+2*i-1] := Complex(0.5);
    end for;
    Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.print(X, 6, "X");
    Modelica_LinearSystems2.Math.Matrices.printMatrix(M,6,"M");
    MM := Complex(1)*M;
    XX := matMul(MM,X);
    XX := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_solve2(
                                                              X, XX);
    Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.print(XX, 6, "XX");
    M:=M-A;
    Modelica_LinearSystems2.Math.Matrices.printMatrix(M,6,"M-A");

  end if;

  K := -Z*transpose(U0)*M;
//  evX := X;
//K:=Modelica_LinearSystems2.StateSpace.Internal.calcK(A,U0,Z,gammaSorted,X);

end assignPolesMI_rob2;
