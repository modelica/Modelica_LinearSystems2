within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zgerq2
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica_LinearSystems2.Math.Complex A[:,:]
    "Square or rectangular matrix";
  output Complex RQ[size(A, 1),size(A, 2)] "RQ factorization in packed format";
  output Complex tau[:] "elementory reflectors";

protected
  Integer n=size(A, 1);
  Integer m=size(A, 2);
  Integer info;
  Real A_real[size(A, 1),size(A, 2)]=A[:, :].re "Square or rectangular matrix";
  Real A_imag[size(A, 1),size(A, 2)]=A[:, :].im "Square or rectangular matrix";
  Real RQ_real[size(A, 1),size(A, 2)]
    "RQ factorization in packed format, real part";
  Real RQ_imag[size(A, 1),size(A, 2)]
    "RQ factorization in packed format, imag part";
  Real tau_real[min(size(A, 1), size(A, 2))]
    "The scalar factors of the elementary reflectors of Q, real part";
  Real tau_imag[min(size(A, 1), size(A, 2))]
    "The scalar factors of the elementary reflectors of Q, imaginary part";

algorithm
  (RQ_real,RQ_imag,tau_real,tau_imag,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zgerq2(
                                                                                            A_real, A_imag);
  for l1 in 1:n loop
    for l2 in 1:m loop
      RQ[l1, l2] := Complex(RQ_real[l1, l2], RQ_imag[l1, l2]);
    end for;
    end for;

   tau := fill(Complex(0), size(tau_real,1));
   for i in 1:min(m, n) loop
     tau[i] := Complex(tau_real[i], tau_imag[i]);
   end for;
end zgerq2;
