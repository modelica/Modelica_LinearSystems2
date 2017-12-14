within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function conversioToTransferFunction_bench
  "Example to compute a zeros-and-poles representation of a MIMO system from state space representation"
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;

  input Integer n=20;
  input Integer ll=7;

  output TransferFunction tf1;
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

  TransferFunction tf[1,1];

algorithm
  for i in 1:n loop
    ss.A[i, i] := n - i + 1;
  end for;
  ss.C[1, 1:ll] := fill(1, ll);

  tf := StateSpace.Conversion.toTransferFunction(ss, 1e-4);
  for i1 in 1:size(ss.C, 1) loop
    for i2 in 1:size(ss.B, 2) loop
      Modelica.Utilities.Streams.print("TransferFunction[" + String(i1) + ","
         + String(i2) + "] = " + String(tf[i1, i2], 12));
    end for;
  end for;

  tf1 := tf[1, 1];
  ok := true;
  annotation(__Dymola_interactive=true);
end conversioToTransferFunction_bench;
