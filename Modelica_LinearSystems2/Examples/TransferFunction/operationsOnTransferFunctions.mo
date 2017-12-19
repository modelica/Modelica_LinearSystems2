within Modelica_LinearSystems2.Examples.TransferFunction;
function operationsOnTransferFunctions
  "Example demonstrating the usage of the functions of Modelica_LinearSystems2.TransferFunction"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Complex;

  output Boolean ok;

protected
  Complex j = Modelica.ComplexMath.j;
  TransferFunction tf1=TransferFunction(n={1,2}, d={2,3,4});
  TransferFunction tf2=TransferFunction(3.4);
  ZerosAndPoles zp1=ZerosAndPoles({-1+0*j},{1+0*j,2+3*j,2-3*j},  k=4);

  ZerosAndPoles zp2=ZerosAndPoles(
    fill(Complex(0), 0),{0.1+0*j}, k=5);
  TransferFunction tf4=ZerosAndPoles.Conversion.toTransferFunction(zp1);
  TransferFunction tf5=ZerosAndPoles.Conversion.toTransferFunction(zp2);

  TransferFunction tf5a=TransferFunction({1}, {1,1});
  TransferFunction tf6=Modelica_LinearSystems2.TransferFunction.s();
  TransferFunction tf7=-tf1;
  TransferFunction tf8=tf1 + tf4;

  TransferFunction tf9=tf1 - tf4;

  TransferFunction tf10=tf1*tf4;
  TransferFunction tf11=tf1/tf4;
  TransferFunction tf12=tf1^2;
  Boolean same1=tf1 == tf1;
  Boolean same2=tf1 == tf4;
  Complex y1=TransferFunction.Analysis.evaluate(tf1, 2+3*j);
  Complex numZeros[:];
  Complex denZeros[:];
  Real k;
  Real k5;
algorithm
  (numZeros,denZeros,k) := TransferFunction.Analysis.zerosAndPoles(tf4);
  print("      tf1 = " + String(tf1));
  print("      tf2 = " + String(tf2));
  print("      tf4 = " + String(tf4));
  (numZeros,denZeros,k5) := TransferFunction.Analysis.zerosAndPoles(tf5);
  tf5a := ZerosAndPoles.Conversion.toTransferFunction(ZerosAndPoles(
    numZeros,
    denZeros,
    k));
  print("      tf5 = " + String(tf5a));
  print("      tf6 = " + String(tf6));
  print("     -tf1 = " + String(tf7));
  print("  tf1+tf4 = " + String(tf8));
  print("  tf1-tf4 = " + String(tf9));
  print("  tf1*tf4 = " + String(tf10));
  print("  tf1/tf4 = " + String(tf11));
  print("  tf1^2   = " + String(tf12));
  print(" tf1==tf1 = " + String(same1));
  print(" tf1==tf4 = " + String(same2));
  print("tf1(2+3j) = " + String(y1));
  print("        k = " + String(k) + " // (zeros,poles,k) = zerosAndPoles(tf4)");
  print("numerator   degree of tf1 = " + String(
    TransferFunction.Analysis.numeratorDegree(tf1)));
  print("denominator degree of tf1 = " + String(
    TransferFunction.Analysis.denominatorDegree(tf1)));
  ok := true;
  annotation (Documentation(info="<html>
<p>
This example shows how to apply various functions from a 
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction\">transfer function</a> portfolio.
</p>
</html>"));
end operationsOnTransferFunctions;
