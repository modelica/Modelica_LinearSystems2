within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function conversionToZerosAndPoles_bench
  "Example to compute a zeros-and-poles representation of a MIMO system from state space representation"
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;

  input Integer n=20;
  input Integer ll=7;

  output ZerosAndPoles zp1;
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

  Modelica_LinearSystems2.ZerosAndPoles zp;

algorithm
  for i in 1:n loop
    ss.A[i, i] := n - i + 1;
  end for;
  ss.C[1, 1:ll] := fill(1, ll);

  zp := StateSpace.Conversion.toZerosAndPoles(ss);
      Modelica.Utilities.Streams.print("ZerosAndPoles-TransferFunction = " + String(zp, 18));

  zp1 := zp;
//  StateSpace.Analysis.analysis(ss);

  ok := true;

end conversionToZerosAndPoles_bench;
