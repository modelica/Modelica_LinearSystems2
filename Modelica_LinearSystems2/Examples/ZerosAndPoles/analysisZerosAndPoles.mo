within Modelica_LinearSystems2.Examples.ZerosAndPoles;
encapsulated function analysisZerosAndPoles "Compute zeros and poles of a ZerosAndPoles transfer function"
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
  Complex numeratorZeros2[4]={-1+j,-1-j,1+0*j,1+0*j};
  Complex denominatorZeros2[6]={1+0*j,2+0*j,2+3*j,2-3*j,3+4*j,3-4*j};

  ZerosAndPoles zp1 = ZerosAndPoles(
    z=numeratorZeros1,
    p=denominatorZeros1);
  ZerosAndPoles zp2 = ZerosAndPoles(
    z=numeratorZeros2,
    p=denominatorZeros2);
  Complex zeros1[:] "Vector of zeros";
  Complex poles1[:] "Vector of poles";
  Real k;

algorithm
  (zeros1,poles1,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp1);

  print("ZerosAndPoles transfer function zp1 = " + String(zp1));
  print("ZerosAndPoles transfer function zp2 = " + String(zp2));
  printComplexVector("numeratorZeros1", numeratorZeros1);
  printComplexVector("\nzeros1", zeros1);

  printComplexVector("\ndenominatorZeros1", denominatorZeros1);
  printComplexVector("\npoles1", poles1);

  ok := true;
end analysisZerosAndPoles;
