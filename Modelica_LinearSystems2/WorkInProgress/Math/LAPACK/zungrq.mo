within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zungrq
  import Complex;

  input Complex RQ[:,:]
    "Square or rectangular matrix";
  input Complex tau[:] "elementary reflectors";
  output Complex Q[size(RQ, 1),size(RQ, 2)] "matrix Q";
protected
  Integer m=size(RQ, 1);
  Integer n=size(RQ, 2);

  Real RQ_real[size(RQ, 1),size(RQ, 2)]=RQ[:,:].re
    "RQ factorization in packed format, real part";
  Real RQ_imag[size(RQ, 1),size(RQ, 2)]=RQ[:,:].im
    "RQ factorization in packed format, imaginary part";
  Real tau_real[min(size(RQ, 1), size(RQ, 2))]=tau[:].re
    "The scalar factors of the elementary reflectors of Q, real part";
  Real tau_imag[min(size(RQ, 1), size(RQ, 2))]=tau[:].im
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  Real Q_real[size(RQ, 1),size(RQ, 2)]
    "matrix Q from RQ factorization in packed format, real part";
  Real Q_imag[size(RQ, 1),size(RQ, 2)]
    "matrix Q from RQ factorization in packed format, imaginary part";

algorithm
  assert(m <= n, "\nInput matrix A[" + String(m) + "," + String(n)
     +"] has more rows as columns. This is not allowed when calling Modelica.Matrices.C_RQ(A).");

  (Q_real,Q_imag) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zungrq(
                                                                                 RQ_real, RQ_imag, tau_real, tau_imag);
  for l1 in 1:m loop
    for l2 in 1:n loop
      Q[l1, l2] := Complex(Q_real[l1, l2],Q_imag[l1, l2]);
    end for;
  end for;

end zungrq;
