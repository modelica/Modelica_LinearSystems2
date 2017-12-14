within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function multiplyRealMatrices3 "Some case studies to complex matrices"
  import Modelica_LinearSystems2;

  input Integer n=50;
  input Integer nr=50;
  output Boolean ok;

protected
  Real C1[n,n]=fill(1,n,n);
  Real C2[n,n];

algorithm
  C2 := Modelica_LinearSystems2.WorkInProgress.Tests.Internal.multiplyRealMatrices3cCall(
                                n,nr,C1);
//  Modelica_LinearSystems2.Math.Matrices.printMatrix(C2,6,"C2");
  ok := true;

end multiplyRealMatrices3;
