within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function realMatrices "Some case studies to real matrices"
  import Modelica_LinearSystems2;

  input Integer n=200;
  output Boolean ok;
protected
  Real C1[n,n]=fill(1,n,n);
  Real C2[n,n];

algorithm
  C2 := C1;
  ok := true;
end realMatrices;
