within Modelica_LinearSystems2.Examples.StateSpace;
function analysisInvariantZeros "Compute invariant zeros of transfer function"
  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input Complex z[:]={-2 + 0*j,-3 + 4*j,-3 - 4*j}
    "Zeros (Complex vector of numerator zeros)";
  input Complex p[:]={-0.5 + 0*j,-5 + 2*j,-5 - 2*j}
    "Poles (Complex vector of denominator zeros)";
  input Real k=1.0 "Constant multiplied with transfer function";

protected
  input Complex j=Modelica_LinearSystems2.Math.Complex.j();

  ZerosAndPoles zp=ZerosAndPoles(
      z=z,
      p=p,
      k=k);

  StateSpace ss=StateSpace(zp);
  Complex Zeros[:];
  Boolean ok;
algorithm
  Zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);

  if size(Zeros, 1) == 0 then
    print("\nSystem\n  "+String(zp)+"\nhas no invariant zeros\n");
  else
    print("\nSystem\n  zp = "+String(zp)+"\n has " + String(size(Zeros, 1)) + " invariant zeros:");
    for i in 1:size(Zeros, 1) loop
      print("   " + String(Zeros[i]));
    end for;
  end if;
  ok := true;
end analysisInvariantZeros;
