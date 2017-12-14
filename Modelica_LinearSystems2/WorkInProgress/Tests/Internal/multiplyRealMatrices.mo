within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function multiplyRealMatrices "Some case studies to complex matrices"
  import Modelica_LinearSystems2;

  input Integer n=30;
  output Boolean ok;
protected
  Real C1[n,n]=fill(1,n,n);
  Real C2[n,n];
  Real t1;
  Real t2;

algorithm
  C2 := C1*C1;
//  Modelica_LinearSystems2.Math.Matrices.printMatrix(C2,6,"C2");
  ok := true;
end multiplyRealMatrices;
