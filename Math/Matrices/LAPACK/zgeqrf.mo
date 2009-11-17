within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function zgeqrf
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica_LinearSystems2.Math.Complex A[:,:]
    "Square or rectangular matrix";
  output Complex Q[size(A, 1),size(A, 2)] "Square upper triangular matrix";
  output Complex R[size(A, 2),size(A, 2)] "Square upper triangular matrix";

protected
  Integer m=size(A, 1);
  Integer n=size(A, 2);
  Integer info;
  Complex QR[size(A, 1),size(A, 2)] "QR factorization in packed format";
  Real A_real[size(A, 1),size(A, 2)]=A[:, :].re "Square or rectangular matrix";
  Real A_imag[size(A, 1),size(A, 2)]=A[:, :].im "Square or rectangular matrix";
  Real QR_real[size(A, 1),size(A, 2)]
    "QR factorization in packed format, real part";
  Real QR_imag[size(A, 1),size(A, 2)]
    "QR factorization in packed format, imaginary part";
  Real Q_real[size(A, 1),size(A, 2)]
    "QR factorization in packed format, real part";
  Real Q_imag[size(A, 1),size(A, 2)]
    "QR factorization in packed format, imaginary part";
  Real tau_real[min(size(A, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  Real tau_imag[min(size(A, 1), size(A_real, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";
//  Complex tau[:] "The scalar factors of the elementary reflectors of Q";

algorithm
assert(m >= n, "\nInput matrix A[" + String(m) + "," + String(n) + "] has more columns as rows.
This is not allowed when calling Modelica.Matrices.C_QR(A).");
  (QR_real,QR_imag,tau_real,tau_imag,info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.wrapper_zgeqrf(
                                                           A_real, A_imag);
  for l1 in 1:m loop
    for l2 in 1:n loop
      QR[l1, l2] := Complex(QR_real[l1, l2], QR_imag[l1, l2]);
    end for;
  end for;

//  tau := fill(Complex(0), size(tau_real, 1));
//   for i in 1:min(m, n) loop
//     tau[i] := Complex(tau_real[i], tau_imag[i]);
//   end for;

// determine R
  R := fill(
    Complex(0),
    n,
    n);
  for i in 1:min(m, n) loop
    for j in i:n loop
      R[i, j] := QR[i, j];
    end for;
  end for;

  (Q_real,Q_imag, info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.wrapper_zungqr(
                                        QR_real, QR_imag, tau_real, tau_imag);
  for l1 in 1:m loop
    for l2 in 1:n loop
      Q[l1, l2] := Complex(Q_real[l1, l2], Q_imag[l1, l2]);
    end for;
  end for;

end zgeqrf;
