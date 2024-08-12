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
  Complex zeros[:] "Vector of zeros";
  Complex poles[:] "Vector of poles";
  Real k;

algorithm
  (zeros,poles,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp1);
  print("Zeros and poles transfer function zp1 = " + String(zp1));
  printComplexVector("\nnumerator zeros for zp1", numeratorZeros1);
  printComplexVector("zeros of zp1", zeros);
  printComplexVector("\ndenominator zeros for zp1", denominatorZeros1);
  printComplexVector("poles of zp1", poles);

  (zeros,poles,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp2);
  print("\nZeros and poles transfer function zp2 = " + String(zp2));
  printComplexVector("\nnumerator zeros for zp2", numeratorZeros2);
  printComplexVector("zeros of zp2", zeros);
  printComplexVector("\ndenominator zeros for zp2", denominatorZeros2);
  printComplexVector("poles of zp2", poles);

  ok := true;
  annotation (
    Documentation(info="<html>
<p>
Construct two zeros and poles transfer functions from complex vectors of numerator and
denominator zeros. Consequently, compute the zeros and poles of the functions and compare
them with the creating complex vectors.
</p>
</html>"));
end analysisZerosAndPoles;
