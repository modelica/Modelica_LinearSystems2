within Modelica_LinearSystems2.Examples.ZerosAndPoles;
encapsulated function analysisZerosAndPoles
  "Compute zeros and poles of a ZerosAndPoles transfer function"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica.ComplexMath.j;
  import Complex;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Internal.printComplexVector;

  output Boolean ok;
protected
  Complex numeratorZeros1[1]={-1+0*j};
  Complex denominatorZeros1[3]={1+0*j,2+3*j,2-3*j};
  Complex numeratorZeros3[4]={-1+j,-1-j,1+0*j,1+0*j};
  Complex denominatorZeros3[6]={1+0*j,2+0*j,2+3*j,2-3*j,3+4*j,3-4*j};

  ZerosAndPoles zp = ZerosAndPoles(
    z=numeratorZeros1,
    p=denominatorZeros1);
  ZerosAndPoles zp2 = ZerosAndPoles(
    z=numeratorZeros3,
    p=denominatorZeros3);
  Complex numeratorZeros2[:];
  Complex denominatorZeros2[:];
  Real k;

algorithm
  (numeratorZeros2,denominatorZeros2,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp);

  print("ZerosAndPoles-TransferFunction1 = " +
    String(zp));
  print("ZerosAndPoles-TransferFunction2 = " +
    String(zp2));
  printComplexVector("numeratorZeros1", numeratorZeros1);
  printComplexVector("\nnumeratorZeros2", numeratorZeros2);

  printComplexVector("\ndenominatorZeros1", denominatorZeros1);
  printComplexVector("\ndenominatorZeros2", denominatorZeros2);

  ok := true;
end analysisZerosAndPoles;
