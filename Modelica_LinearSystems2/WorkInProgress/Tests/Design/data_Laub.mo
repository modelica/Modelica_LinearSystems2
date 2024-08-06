within Modelica_LinearSystems2.WorkInProgress.Tests.Design;
function data_Laub "Example for pole assignment"
  extends Modelica.Icons.Function;

  import Complex;
  import Modelica_LinearSystems2.WorkInProgress.Tests.Internal.DesignData;

  input Integer n=10 annotation(Dialog);
  input Integer m=1;
  input Real alpha=0.1;
  output DesignData data(
  redeclare Real A[n,n],
  redeclare Real B[n,m],
  redeclare Complex assignedPoles[n]);

protected
  Real A[n,n]=diagonal(array(-(n-i) for i in 1:n));
  Real B[n,m];
  Complex pp[:]=array(Complex(-i-10) for i in 2:2:2*n);

algorithm
  for i in 1:n-1 loop
    A[i+1,i] := alpha;
  end for;
  for i in 1:m loop
    B[i,i] := 1;
  end for;

  data.A:=A;
  data.B:=B;
  data.assignedPoles:=pp;
  data.K:=fill(0, 0, 0);
end data_Laub;
