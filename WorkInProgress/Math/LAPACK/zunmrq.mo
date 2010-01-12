within Modelica_LinearSystems2.WorkInProgress.Math.LAPACK;
function zunmrq
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica_LinearSystems2.Math.Complex C[:,:] "Matrix to multiply with Q";
  input Modelica_LinearSystems2.Math.Complex RQ[:,:]
    "Matrix as a result of zgerq2";
  input Complex tau[:] "elementary reflectors";
  input Boolean left=true;
  input Boolean trans=false;
  output Complex CQ[size(C, 1),size(C, 2)]
    "Matrix product C*Q, Q*C, C*Q**H, Q**H*C";
protected
  Integer m=size(C, 1);
  Integer n=size(C, 2);

  Real C_real[size(C, 1),size(C, 2)]=C[:,:].re
    "Real part of matrix to multiply with Q";
  Real C_imag[size(C, 1),size(C, 2)]=C[:,:].im
    "Imag part of matrix to multiply with Q";
  Real RQ_real[size(RQ, 1),size(RQ, 2)]=RQ[:,:].re
    "RQ factorization in packed format, real part";
  Real RQ_imag[size(RQ, 1),size(RQ, 2)]=RQ[:,:].im
    "RQ factorization in packed format, imaginary part";
  Real tau_real[min(size(RQ, 1), size(RQ, 2))]=tau[:].re
    "The scalar factors of the elementary reflectors of Q, real part";
  Real tau_imag[min(size(RQ, 1), size(RQ, 2))]=tau[:].im
    "The scalar factors of the elementary reflectors of Q, imaginary part";
  Real CQ_real[size(CQ, 1),size(CQ, 2)] "Matrix product, real part";
  Real CQ_imag[size(CQ, 1),size(CQ, 2)] "Matrix product, imaginary part";

algorithm
  (CQ_real,CQ_imag) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.wrapper_zunmrq(
                                                                                   RQ_real, RQ_imag, tau_real, tau_imag, C_real, C_imag, left, trans);
  for l1 in 1:m loop
    for l2 in 1:n loop
      CQ[l1, l2] := Complex(CQ_real[l1, l2],CQ_imag[l1, l2]);
    end for;
  end for;

end zunmrq;
