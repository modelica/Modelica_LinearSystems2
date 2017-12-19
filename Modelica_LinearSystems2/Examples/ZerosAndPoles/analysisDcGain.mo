within Modelica_LinearSystems2.Examples.ZerosAndPoles;
encapsulated function analysisDcGain
  "Example for computation of steady state gain of a zeros and poles representation"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Utilities.Types;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Complex;

  output Boolean ok;
protected
  Complex j = Modelica.ComplexMath.j;
  Complex numeratorZeros1[1]={-1+0*j};
  Complex denominatorZeros1[3]={1+0*j,2+3*j,2-3*j};
  Complex numeratorZeros3[4]={-1+j,-1-j,1+0*j,1+0*j};
  Complex denominatorZeros3[6]={1+0*j,2+0*j,2+3*j,2-3*j,3+4*j,3-4*j};

  ZerosAndPoles zp=ZerosAndPoles(
    z=numeratorZeros1,
    p=denominatorZeros1);
  ZerosAndPoles zp2=ZerosAndPoles(
    z=numeratorZeros3,
    p=denominatorZeros3);
  ZerosAndPoles zp3=ZerosAndPoles.Internal.baseFilter(
    Modelica_LinearSystems2.Utilities.Types.AnalogFilter.Bessel,
    order=5);
  Real k;
  Boolean finite;
algorithm
  ok := false;

  (k, finite) := ZerosAndPoles.Analysis.dcGain(zp);
  print("k1 = " + String(k) + ", finite1 = " + String(finite));

  (k, finite) := ZerosAndPoles.Analysis.dcGain(zp2);
  print("k2 = " + String(k) + ", finite2 = " + String(finite));

  (k, finite) := ZerosAndPoles.Analysis.dcGain(zp3);
  print("k3 = " + String(k) + ", finite3 = " + String(finite));

  ok := true;
  annotation (Documentation(info="<html>
<p>
This example shows how to calculate the
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.dcGain\">steady state gain <b>K</b></a>
of a zeros and poles representation given internally.
</p>
</html>
"));
end analysisDcGain;
