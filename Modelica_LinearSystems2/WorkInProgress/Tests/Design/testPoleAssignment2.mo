within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function testPoleAssignment2
  "Function to assess algorithms for pole assignment"
  extends Modelica.Icons.Function;

  import Complex;
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;
  import Modelica.Utilities.Streams;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Design;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;
  import Modelica_LinearSystems2.StateSpace;

  input String dataFile=TestDataDir + "data_Byers6.mat" annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                     caption="state space system data file"),enable = systemOnFile));
  input Types.AssignPolesMethod method=Tests.Types.AssignPolesMethod.KNV
    "method for pole assignment";
  input Boolean isSI=true;
  input String outputFile = "";
  input Boolean deleteExistingOutputfile=true;

protected
  Integer nm[2]=Streams.readMatrixSize(dataFile, "B")
    "Read system order and number of inputs";
  Integer nmk[2]=Streams.readMatrixSize(dataFile, "K") "Read dimensions of K";
  Real A[:,:]=Streams.readRealMatrix(dataFile, "A", nm[1], nm[1]) "Read system matrix A";
  Real B[:,:]=Streams.readRealMatrix(dataFile, "B", nm[1], nm[2]) "Read system matrix B";
  Complex j = Modelica.ComplexMath.j;
  Real assignedPolesR[1,:]=Streams.readRealMatrix(dataFile, "assignedPoles", 1, nm[1])
    "Read real part of assigned poles";
  Real assignedPolesI[1,:]=Streams.readRealMatrix(dataFile, "assignedPolesIm", 1, nm[1])
    "Read imaginary part of assigned poles";
  Complex assignedPoles[:]=Complex(1)*assignedPolesR[1, :] + j*assignedPolesI[1, :]
    "Complex assigned poles";

  Boolean isKprovided = min(nmk) > 0;
  Real Ki[:,:] = if isKprovided then
    Streams.readRealMatrix(dataFile, "K", nm[2], nm[1]) else fill(0, 0, 0);
//  Integer n=size(A, 1);
  Real S[nm[1],nm[1]] "closed loop system matrix A-BK";
  StateSpace ss = Modelica_LinearSystems2.StateSpace(
    A=A, B=B, C=zeros(1, nm[1]), D=zeros(1, size(B, 2)));
  Real Xre[nm[1],nm[1]];
  Real k[size(A, 1)] "Feedback gain matrix";

public
  output Real K[size(B, 2),size(A, 1)] "Feedback gain matrix";
  output Complex calcPoles[:];
  output Real kappa2 "condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2";
  output Real kappaF "condition number kappa_F(X) = ||X||_F * ||inv(X)||_F";
  output Real zeta
    "condition number by Byers, zeta(X) = (||X||_F)^2 + (||inv(X)||_F)^2";
  output Real cInf "condition number vu1=||c||_inf = max(c_j)";
  output Real nu2 "Euclidean norm of the feedback matrix";
  output Real nuF "Frobenius norm of the feedback matrix";
  output Real dlambda "Distance between the assigned and the calculated poles";
  output Real gap=0.0;
  output Real Jalpha[11]
    "Combined condition number, JKX=alpha/2*(kappa2X_B) + (1-alpha)/2*normFroK^2";
  output Complex X[:,:] "right eigenvectors of the closed loop system";

algorithm
// use single input algorithm
  if isSI and nm[2] == 1 then
    K := Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal.assignPolesSI_rq(ss, assignedPoles);
//    k := Modelica_LinearSystems2.StateSpace.Design.assignPolesSI(ss, assignedPoles);
//    K := transpose(matrix(k));
    ss.A := ss.A - ss.B*K;
    S := ss.A;
    (X,calcPoles) :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenVectors(S);
//    X := Complex(1)*Xre;
  else// end isSI

    if method == Tests.Types.AssignPolesMethod.KNV then
// extented robust KNV-algortihm according to MATLAB's place-function
      (K,X) := Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal.assignPolesMI_rob(
                                                                             A, B, assignedPoles);
      S := A - B*K;
      calcPoles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(S);
      if isKprovided then
        gap := Modelica.Math.Matrices.norm(K - Ki);
      end if;
    elseif method ==Modelica_LinearSystems2.WorkInProgress.Tests.Types.AssignPolesMethod.Schur then
// Schur method
      (K,S,calcPoles,,,,X) :=
        Modelica_LinearSystems2.WorkInProgress.StateSpace.Design.assignPolesMI(ss, assignedPoles, -1e10, Modelica.Math.Matrices.norm(ss.A, 1)*1e-12, true);
//        Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, assignedPoles, -1e10, Modelica.Math.Matrices.norm(ss.A, 1)*1e-12, true);

      if isKprovided then
        gap := Modelica.Math.Matrices.norm(K - Ki);
      end if;
    else
      assert(false, "Argument method (= " + String(method) + ") of testPoleAssignment is wrong. It has to be \"KNV\" or \"Schur\"");
    end if;
  end if;
  // calculate condition numbers
  (kappa2,kappaF,,cInf,nu2,nuF,zeta,Jalpha,dlambda) := conditionNumbers(K, X, assignedPoles, calcPoles);
  if deleteExistingOutputfile then
    if Modelica.Utilities.Files.exist(outputFile) then
      Modelica.Utilities.Files.removeFile(outputFile);
    end if;
    end if;
  if isSI and nm[2] == 1 then
    print("---- Using SI-algorithm ----\n",outputFile);
  else
    print("---- Using MI-algorithm ----\n",outputFile);
  end if;
  print("n = "+String(nm[1])+",  m = "+ String(nm[2])+"\n",outputFile);
  print(Matrices.printMatrix(K, 6, "K"),outputFile);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print(
    "assignedPoles",
    assignedPoles,
    outputFile);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print(
    "calcPoles",
    calcPoles,
    outputFile);
//   Matrices.printMatrix(Re(X), 6, "ReX");
//   Matrices.printMatrix(Im(X), 6, "ImX");
  print("kappa2 " + String(kappa2),outputFile);
  print("kappaF " + String(kappaF),outputFile);
  print("zeta " + String(zeta),outputFile);
  print("cInf " + String(cInf),outputFile);
  print("nu2 " + String(nu2),outputFile);
  print("nuF " + String(nuF),outputFile);
  print("dlambda " + String(dlambda),outputFile);
  if isKprovided then
    print("gap " + String(gap),outputFile);
  end if;
  print("Jalpha = " + Modelica_LinearSystems2.Math.Vectors.printVector(Jalpha),outputFile);
  print("\nMLS2 & $"+ String(kappa2)+"$ & $"+ String(zeta)+"$ & $"+ String(cInf)+"$ & $"+ String(nu2)+"$ & $"+ String(dlambda)+"$\\   \\hline",outputFile);

  annotation (Documentation(info="<html>
<p>
Computes the feedback gain K for the state space system according to assigned close loop poles.
</p>
</html>"));
end testPoleAssignment2;
