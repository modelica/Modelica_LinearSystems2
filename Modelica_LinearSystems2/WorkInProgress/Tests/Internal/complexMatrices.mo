within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function complexMatrices "Some case studies to complex matrices"
  import Complex;

  input Integer n=200;
  output Boolean ok;
protected
  Complex C1[n,n]=fill(Complex(1,1),n,n);
  Complex C2[n,n];

algorithm
  C2 := C1;
  ok := true;
end complexMatrices;
