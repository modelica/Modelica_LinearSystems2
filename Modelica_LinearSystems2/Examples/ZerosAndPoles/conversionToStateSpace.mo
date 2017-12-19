within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function conversionToStateSpace
  "Example for transformation of a transfer function from zeros and poles representation into a StateSpace description"
  extends Modelica.Icons.Function;

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

  annotation (Documentation(info="<html>
<p>
This example shows how to transform a transfer function from zeros and poles representation
<code>zpi</code> into a state-space representation by either using the explicit conversion
<p>
<blockquote><pre>
ZerosAndPoles.Conversion.toStateSpace(zpi);
</pre></blockquote>
<p>
or the short conversion with overloadig
</p>
<blockquote><pre>
StateSpace ss2=StateSpace(zpi);
</pre></blockquote>
</html>"));
end conversionToStateSpace;
