within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zungqr
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica_LinearSystems2.Math.Complex QR[:,:]
    "Square or rectangular matrix";
  input Complex tau[:] "elementary reflectors";
  output Complex Q[size(QR, 1),size(QR, 2)] "matrix Q";
protected
  Integer m=size(QR, 1);
  Integer n=size(QR, 2);

  Real QR_real[size(QR, 1),size(QR, 2)]=QR[:,:].re
    "QR factorization in packed format, real part";
  Real QR_imag[size(QR, 1),size(QR, 2)]=QR[:,:].im
    "QR factorization in packed format, imaginary part";
  Real tau_real[min(size(QR, 1), size(QR, 2))]=tau[:].re
    "The scalar factors of the elementary reflectors of Q, real part";
  Real tau_imag[min(size(QR, 1), size(QR, 2))]=tau[:].im
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  Real Q_real[size(QR, 1),size(QR, 2)]
    "matrix Q from QR factorization in packed format, real part";
  Real Q_imag[size(QR, 1),size(QR, 2)]
    "matrix Q from QR factorization in packed format, imaginary part";

algorithm
  (Q_real,Q_imag) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zungqr(
                                  QR_real, QR_imag, tau_real, tau_imag);
  for l1 in 1:m loop
    for l2 in 1:n loop
      Q[l1, l2] := Complex(Q_real[l1, l2],Q_imag[l1, l2]);
    end for;
  end for;

end zungqr;
