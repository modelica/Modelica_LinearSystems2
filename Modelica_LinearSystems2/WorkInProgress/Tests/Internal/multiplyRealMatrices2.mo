within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function multiplyRealMatrices2 "Some case studies to complex matrices"
  import Modelica_LinearSystems2;

  input Integer n=50;
  input Integer nr=50;
  output Boolean ok;
protected
  Real C1[n,n]=fill(1,n,n);
  Real C2[n,n]=C1;

algorithm
  for i in 1:nr loop
   C2 := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgemm(C1,C2,C1);
  end for;

//  Modelica_LinearSystems2.Math.Matrices.printMatrix(C2,6,"C2");
  ok := true;
end multiplyRealMatrices2;
