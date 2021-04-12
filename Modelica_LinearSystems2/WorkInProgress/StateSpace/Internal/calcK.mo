within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function calcK
  "Computes the feedback matrix from the assigned eigenvalues, closed loop eigenvectors and the B matrix factorization"
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;

  input Real A[:,size(A, 1)] "Real square system matrix";
  input Real U0[size(A, 1),:] "U0 and Z are the decompositions of B";
  input Real Z[size(U0, 2),size(U0, 2)] "Z and U0 are the decompositions of B";
  input Complex gamma[size(A, 1)] "Assigned complex eigenvalues";
  input Complex X[size(A, 1),size(A, 1)] "Closed loop eigenvectors";
  input Integer nre "number of real eigenvalues";

  output Real K[size(U0, 2),size(A, 1)] "Feedback matrix";

protected
  Integer n=size(A, 1);
  Integer m=size(U0, 2);

  Real gamma_real[n]=Re(gamma) "Eigenvalue vector, real part";
  Real gamma_imag[n]=Im(gamma) "Eigenvalue vector, imaginary part";
  Real X_real[n,n]=Re(X) "Eigenvectors, real part";
  Real X_imag[n,n]=Im(X) "Eigenvectors, imaginary part";

algorithm
  K := Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal.wrapper_calcK(
                                                                 A, U0, Z, gamma_real, gamma_imag, X_real,X_imag, nre);
end calcK;
