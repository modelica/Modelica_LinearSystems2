within Modelica_LinearSystems2.WorkInProgress.Tests.Complex;
model ComplexOperations
  import Complex;
  import Modelica.ComplexMath.j;
  Complex c1;
  Complex c3;
  Complex c4;
equation
  c1 = Complex(2);
  c3 = 2+3*j;
  c4 = c3 + c1*j + j/c1;
end ComplexOperations;
