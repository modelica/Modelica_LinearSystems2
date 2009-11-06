within Modelica_LinearSystems2.Tests.Design;
function poleAssignmentKautsky2 "Example for pole assignment"
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

  Real A[5,5]=[-0.1094,  0.0628, 0.0,     0.0,  0.0;
                1.306,  -2.132,   0.9807,  0.0,  0.0;
                0.0,  1.595,   -3.149,  1.547,  0.0;
                0.0,  0.0355,  2.632,  -4.257, 1.855;
                0.0,  0.00227,  0.0,  0.1636,  -0.1625];

  Real B[5,2]=[0.0,  0.0;  0.0638,  0.0;  0.0838,  -0.1396;  0.1004,  -0.206;  0.0063,  -0.0128];

  Real C[1,5]=[1,0,0,0,0];
  Real D[1,2]=[0,0];
  Integer n=size(A, 1);

  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=A[1:n, 1:n],
      B=B[1:n, :],
      C=C[:, 1:n],
      D=D);

  Modelica_LinearSystems2.Math.Complex assignedPoles[:]=cat(1,Complex(1)*{-0.2,-0.5,-1.0},{-1.0+j, -1.0-j});
//  Modelica_LinearSystems2.Math.Complex assignedPoles[:]=cat(1,Complex(1)*{-0.2},{-2+j,-2-j,-1.0+j, -1.0-j});
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
  // extented robust KNV-algortihm according to MATLAB's place-function
  (K, X) := Modelica_LinearSystems2.StateSpace.Internal.assignPolesMI_rob5(ss.A, ss.B, assignedPoles);
  Matrices.printMatrix(K,6,"K");
  S := ss.A - ss.B*K;
  newPoles := Complex.eigenValues(S);
  Modelica_LinearSystems2.Math.Complex.Vectors.print("assignedPoles", assignedPoles);
  Complex.Vectors.print("newPoles", newPoles);
  (kappa2X, kappaFroX, kappaFroYT, cInf, norm2K, normFroK, kappa2X_B, JXK) :=  conditionNumbers(K, X, assignedPoles, newPoles);

end poleAssignmentKautsky2;
