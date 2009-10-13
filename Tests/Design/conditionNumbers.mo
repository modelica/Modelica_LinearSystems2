within Modelica_LinearSystems2.Tests.Design;
function conditionNumbers
  "Calculate several condition numbers to evaluate a pole assigment method"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input Real K[:,:] "state feedback matrix";
  input Complex X[:,:] "right eigenvectors of the closed loop system";
  output Real kappa2X "condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2";
  output Real kappaFroX "condition number kappa_F(X) = ||X||_F * ||inv(X)||_F";
  output Real kappaFroYT
    "condition number kappa_F(YT) = ||YT||_F * ||inv(YT)||_F";
  output Real cInf "condition number vu1=||c||_inf = max(c_j)";
  output Real norm2K "Euclidean norm of the feedback matrix";
  output Real normFroK "Frobenius norm of the feedback matrix";
  output Real kappa2X_B
    "condition number by Byers, kappa_2(X) = (||X||_2)^2 + (||inv(X)||_2)^2";
  output Real JXK[11]
    "condition number by Varga, JKX=alpha/2*(kappa2X_B) + (1-alpha)/2*normFroK^2";
protected
  Integer n=size(X,1);
  Real c[size(X, 1)] "condition number cj";
  Complex YT[:,:] "left eigenvectors of the closed loop system";

algorithm
  YT := Complex.Matrices.inv(X);
  for l1 in 1:n loop
    c[l1] := Complex.Vectors.norm(YT[l1, :]);
  end for;
  //performance indices
  // condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2
  kappa2X := Complex.Matrices.conditionNumber(X);
  print("kappa2X = " + String(kappa2X));
  kappaFroX := Complex.Matrices.conditionNumber(X, 3);
  print("kappaF = " + String(kappaFroX));
  cInf := Modelica.Math.Vectors.norm(c, Modelica.Constants.inf);
  print("cInf = " + String(cInf));
  norm2K := Modelica.Math.Matrices.norm(K);
  print("norm2K = " + String(norm2K));
  normFroK := Matrices.norm(K, 3);
  print("normFroK = " + String(normFroK));
  kappaFroYT := Complex.Matrices.conditionNumber(YT, 3);
  kappa2X_B := kappaFroX^2 + kappaFroX^2;
  print("kappa2X_B = " + String(kappa2X_B));
  JXK := array(alpha/2*kappa2X_B + (1 - alpha)/2*normFroK^2 for alpha in 0:0.1:1);
  Math.Vectors.printVector(JXK, 6, "JXK");
end conditionNumbers;
