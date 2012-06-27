within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function conversionToStateSpace
  "Transform a transfer function from zeros and poles representation into a StateSpace description"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input ZerosAndPoles zpi=ZerosAndPoles(k=1, n1={1}, n2=fill(0,0,2),d1=fill(0,0),d2=[1,1]);

  output Boolean ok;

protected
  StateSpace ss1=ZerosAndPoles.Conversion.toStateSpace(zpi);      //explicit conversion
  StateSpace ss2=StateSpace(zpi);                                 //short conversion using overloadig
  ZerosAndPoles zpo1=StateSpace.Conversion.toZerosAndPoles(ss1);  //explicit conversion

algorithm
   Modelica.Utilities.Streams.print("zpi = " + String(zpi));
   Modelica.Utilities.Streams.print("ss1 = " + String(ss1,6,"ss1"));
   Modelica.Utilities.Streams.print("ss2 = " + String(ss2,6,"ss2"));
   Modelica.Utilities.Streams.print("zpo1 = " + String(zpo1));

   ok := true;

end conversionToStateSpace;
