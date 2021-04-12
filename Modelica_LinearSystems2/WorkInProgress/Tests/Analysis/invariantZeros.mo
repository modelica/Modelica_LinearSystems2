within Modelica_LinearSystems2.WorkInProgress.Tests.Analysis;
function invariantZeros
  "Example to compute the invariant zeros of a state space system"
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;
  import Complex;

  input Integer n=20;
  input Integer ll=7;

  output Complex iz[:];
  output Boolean ok;

protected
  Real A[n,n];
  Real B[n,1]=fill(
      1,
      n,
      1);
  Real C[1,n]=fill(
      0,
      1,
      n);
  Real D[1,1]=[0];

  StateSpace ss=StateSpace(
      A=A,
      B=B,
      C=C,
      D=D);

algorithm
  for i in 1:n loop
    ss.A[i, i] := n - i + 1;
  end for;
  ss.C[1, 1:ll] := fill(1, ll);

  iz := StateSpace.Analysis.invariantZeros(ss);
  ok := true;

end invariantZeros;
