within Modelica_LinearSystems2.ComplexMathAdds.Examples;
function addTwoComplexNumbers "Show how to add 2 complex number"
  extends Modelica.Icons.Function;
  import Modelica.Utilities.Streams;
  import Modelica.ComplexMath.j;

  output Boolean ok;
protected
  Complex c1=2+3*j;
  Complex c2=3+4*j;
  Complex c3;
algorithm

  Streams.print("c1 = " + String(c1));
  Streams.print("c2 = " + String(c2));

  c3 := c1 + c2;
  Streams.print("c3 = c1 + c2 = " + String(c3));
  ok := true;
end addTwoComplexNumbers;
