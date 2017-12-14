within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function conditionNumbers
  "Calculate several condition numbers to evaluate a pole assignment method"
  extends Modelica.Icons.Function;

  import Complex;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input Real K[:,:] "state feedback matrix";
  input Complex X[:,:] "right eigenvectors of the closed loop system";
  input Complex assignedPoles[size(X,1)]=fill(Complex(0),size(X,1));
  input Complex calcPoles[size(X,1)]=fill(Complex(0),size(X,1));
  output Real kappa2X "condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2";
  output Real kappaFroX "condition number kappa_F(X) = ||X||_F * ||inv(X)||_F";
  output Real kappaFroYT
    "condition number kappa_F(YT) = ||YT||_F * ||inv(YT)||_F";
  output Real cInf "condition number vu1=||c||_inf = max(c_j)";
  output Real norm2K "Euclidean norm of the feedback matrix";
  output Real normFroK "Frobenius norm of the feedback matrix";
  output Real kappa2X_B
    "condition number by Byers, kappa_2XB(X) = (||X||_F)^2 + (||inv(X)||_F)^2";
  output Real JXK[11]
    "condition number by Varga, JKX=alpha*(kappa2X_B)/2 + (1-alpha)*normFroK^/22";
  output Real dl
    "gap between the assigned and the calculated poles dl=norm(ap-cp)/max(1,norm(ap))";
protected
  Integer n=size(X,1);
  Real c[size(X, 1)] "condition number cj";
  Complex YT[:,:] "left eigenvectors of the closed loop system";
  Complex sortedAssignedPoles[size(X,1)];
  Complex sortedCalcPoles[size(X,1)];

algorithm
  sortedAssignedPoles := Modelica.ComplexMath.Vectors.sort(assignedPoles);
  sortedCalcPoles := Modelica.ComplexMath.Vectors.sort(calcPoles);
  dl := Modelica.ComplexMath.Vectors.norm(sortedAssignedPoles - sortedCalcPoles)/max(1, Modelica.ComplexMath.Vectors.norm(sortedAssignedPoles));
  YT := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.inv(X);
  for l1 in 1:n loop
    c[l1] := Modelica.ComplexMath.Vectors.norm(YT[l1, :])*Modelica.ComplexMath.Vectors.norm(X[:, l1])/Modelica.ComplexMath.'abs'(Complex.'*'.scalarProduct(YT[l1, :], X[:, l1]));
  end for;
  //performance indices
  // condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2
  kappa2X := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.conditionNumber(
                                              X);
  kappaFroX := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.conditionNumber(
                                                X, 3);
  cInf := Modelica.Math.Vectors.norm(c, Modelica.Constants.inf);
  norm2K := Modelica.Math.Matrices.norm(K);
  normFroK := Matrices.norm(K, 3);
  kappaFroYT := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.conditionNumber(
    YT, 3);
  kappa2X_B := Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.norm(
    X)^2 + Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices.norm(
    YT)^2;
  JXK := array(alpha/2*kappa2X_B + (1 - alpha)/2*normFroK^2 for alpha in 0:0.1:1);
  annotation (Documentation(info="<html>
</html>"));
end conditionNumbers;
