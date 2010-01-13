within Modelica_LinearSystems2.WorkInProgress.Tests.Complex;
model ComplexOperations
  import Modelica_LinearSystems2.Math.Complex;
  Complex c1;
  Complex c3;
  Complex c4;
  Complex j;
equation
  c1 = Complex(2);
  j  = Complex.j();
  c3 = 2+3*j;
  c4 = c3 + c1*j + j/c1;
end ComplexOperations;
