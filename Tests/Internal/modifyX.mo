within Modelica_LinearSystems2.Tests.Internal;
function modifyX
  "Contains a C sub routine of robust pole assignment to modify the eigenvector matrix X according to Kautsky algorithm"

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Complex;
  import Re = Modelica_LinearSystems2.Math.Complex.real;
  import Im = Modelica_LinearSystems2.Math.Complex.imag;

  input Complex X[:,size(X,1)] "Complex eigenvector matrix";
  input Complex S[size(X,1),:] "Complex eigenvector matrix";
  input Integer m
    "Rank of the system input matrix B; S_real and S_imag must have n*m columns";
  input Integer ncp "number of complex pairs";
  input Integer steps "Number of iterations";

  output Complex Xm[size(X, 1),size(X, 2)];

protected
   Complex j=Modelica_LinearSystems2.Math.Complex.j();
   Integer n=size(X,1);
   Real X_real[n,n]=Re(X) "Eigenvector matrix, real part";
   Real X_imag[n,n]=Im(X) "Eigenvector matrix, imaginary part";
   Real S_real[n,m*n]=Re(S) "Eigenvector bases, real part";
   Real S_imag[n,m*n]=Im(S) "Eigenvector bases, imaginary part";

  Real Xm_real[n,n];
  Real Xm_imag[n,n];

  Integer i;
  Integer ii;

algorithm
  (Xm_real, Xm_imag) :=Modelica_LinearSystems2.StateSpace.Internal.wrapper_modifyX(X_real, X_imag, n, S_real, S_imag, m, ncp, steps);
  for i in 1:n loop
    for ii in 1:n loop
      Xm[i,ii] := Complex(Xm_real[i,ii],Xm_imag[i,ii]);
    end for;
  end for;
end modifyX;
