within Modelica_LinearSystems2.Tests.Design;
function poleAssignmentByers5 "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Re = Modelica_LinearSystems2.Math.Complex.real;
  import Im = Modelica_LinearSystems2.Math.Complex.imag;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;
  annotation (Documentation(info="<html>
<p>
Computes the gain vector k for the state space system
<pre> 
ss = StateSpace(A=[-1,1;0,-2],B=[0, 1],C=[1,0; 0, 1],D=[0; 0])
</pre>
such that for the state feedback 
<pre>u = -k*y = -k*x</pre> the closed-loop
poles are placed at 
<pre>p = {-3,-4}.</pre>
</html>"));

protected
  Complex j = Modelica_LinearSystems2.Math.Complex.j();
  Modelica_LinearSystems2.Math.Complex newPoles[:];
  Real S[:,:] "closed loop system matrix A-BK";

  Real A[5,5]=diagonal({1, 10, 0.1, 0.1, 10})*[-1.29e-1,  0.0,  3.96e-2, 2.5e-2,   1.91e-2;
                                                3.29e-3,  0.0, -7.79e-5, 1.22e-4, -6.21e-1;
                                                7.18e-2,  0.0, -1.0e-1,  8.87e-4,  -3.85;
                                                4.11e-2,  0.0,  0.0,    -8.22e-2,   0.0;
                                                3.51e-4,  0.0,  3.5e-5,  4.26e-5,  -7.43e-2]*diagonal({1, 0.1, 10, 10, 0.1});

  Real B[5,2]=diagonal({1, 10, 0.1, 0.1, 10})*[0.0,  1.39e-3; 0.0,  3.59e-5; 0.0,  -9.89e-3; 2.49e-5,  0.0; 0.0,  -5.34e-6]*diagonal({10000, 100});

  Real C[1,5]=[1,0,0,0,0];
  Real D[1,2]=[0,0];
  Integer n=size(A, 1);

  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=A[1:n, 1:n],
      B=B[1:n, :],
      C=C[:, 1:n],
      D=D);

  Modelica_LinearSystems2.Math.Complex assignedPoles[:]=Complex(1)*{-0.01,  -0.02,  -0.03,  -0.04,  -0.05};
  Real K[size(B, 2),size(A, 1)] "Feedback gain matrix";

//  Real c[size(A,1)] "condition number cj";
  Real kappa2X "condition number kappa_2(X) = ||X||_2 * ||inv(X)||_2";
  Real kappaFroX "condition number kappa_F(X) = ||X||_F * ||inv(X)||_F";
  Real kappaFroYT "condition number kappa_F(YT) = ||YT||_F * ||inv(YT)||_F";
  Real cInf "condition number vu1=||c||_inf = max(c_j)";
  Real norm2K "Euclidean norm of the feedback matrix";
  Real normFroK "Frobenius norm of the feedback matrix";
  Real kappa2X_B
    "condition number by Byers, kappa_2(X) = (||X||_2)^2 + (||inv(X)||_2)^2";
  Real JXK[11]
    "condition number by Varga, JKX=alpha/2*(kappa2X_B) + (1-alpha)/2*normFroK^2";

  Complex X[:,:] "right eigenvectors of the closed loop system";
//  Complex YT[:,:] "left eigenvectors of the closed loop system";
  Real evec[size(A,1),size(A,1)] "eigenvectors of A-BK";

  Integer l1;

  Real alpha;

algorithm
  // Schur method Varga
  (K,S,newPoles,,,,X) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, assignedPoles);
  Modelica_LinearSystems2.Math.Matrices.printMatrix(K,6,"K");
  Modelica_LinearSystems2.Math.Complex.Vectors.print("assignedPoles", assignedPoles);
  Modelica_LinearSystems2.Math.Complex.Vectors.print("newPoles", newPoles);
  Matrices.printMatrix(Re(X),6,"ReX");
  Matrices.printMatrix(Im(X),6,"ImX");
  (,evec) := Modelica.Math.Matrices.eigenValues(S);
  Matrices.printMatrix(evec,6,"evec");
  (kappa2X, kappaFroX, kappaFroYT, cInf, norm2K, normFroK, kappa2X_B, JXK) :=  conditionNumbers(K, X);

  // extented robust KNV-algortihm according to MATLAB's place-function
  (K, X) := Modelica_LinearSystems2.StateSpace.Internal.assignPolesMI_rob(ss.A, ss.B, assignedPoles);
  Matrices.printMatrix(K,6,"K");
  S := ss.A - ss.B*K;
  newPoles := Complex.eigenValues(S);
  Modelica_LinearSystems2.Math.Complex.Vectors.print("assignedPoles", assignedPoles);
  Complex.Vectors.print("newPoles", newPoles);
  (kappa2X, kappaFroX, kappaFroYT, cInf, norm2K, normFroK, kappa2X_B, JXK) :=  conditionNumbers(K, X);

end poleAssignmentByers5;
