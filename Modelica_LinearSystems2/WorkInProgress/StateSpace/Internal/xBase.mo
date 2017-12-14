within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function xBase "Compute the eigenvector bases according to Kautsky algorithm"
  import Modelica_LinearSystems2;
  import Complex;
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;

  input Real A[:,size(A,1)] "Real square system matrix";
  input Real B[size(A,1),:] "Real input matrix";
  input Complex gamma[size(A,1)] "Assigned complex eigenvalues";
  input Integer ncp "Number of complex pairs of eigenvalues";

  output Real U0[size(A, 1),size(B, 2)] "U0 and Z are the decompositions of B";
  output Real Z[size(B,2),size(B,2)] "Z and U0 are the decompositions of B";
  output Complex S[size(A,1),(size(A,1)-ncp)*size(B,2)] "Eigenvector bases";
  output Integer rankB;

protected
  Complex j = Modelica.ComplexMath.j;
  Integer n=size(A,1);
  Integer m=size(B,2);

  Real gamma_real[n]=Modelica.ComplexMath.real(
                         gamma) "Eigenvalue vector, real part";
  Real gamma_imag[n]=Modelica.ComplexMath.imag(
                         gamma) "Eigenvalue vector, imaginary part";
  Real S_real[n,m*(n-ncp)] "Eigenvector bases, real part";
  Real S_imag[n,m*(n-ncp)] "Eigenvector bases, imaginary part";
  Integer i;
  Integer ii;

algorithm
  (U0, Z, S_real, S_imag, rankB) :=Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal.wrapper_xBase(
                                                                                            A, B, gamma_real, gamma_imag, ncp);
  assert(m==rankB,"Input matrix B must have full column rank");
  for i in 1:n loop
    for ii in 1:(n-ncp)*m loop
      S[i,ii] := Complex(S_real[i,ii],S_imag[i,ii]);
    end for;
  end for;
end xBase;
