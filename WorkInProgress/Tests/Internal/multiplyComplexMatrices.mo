within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function multiplyComplexMatrices "Some case studies to complex matrices"
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2;

  input Integer n=30;
  output Boolean ok;
protected
  Complex C1[n,n]=fill(Complex(1,1),n,n);
  Complex C2[n,n];
  Real t1;
  Real t2;

algorithm
//  C2 := C1*C1;
  C2 := Complex.Matrices.matMatMul(C1,C1);
  ok := true;
end multiplyComplexMatrices;
