within Modelica_LinearSystems2.Examples.StateSpace;
function conversionFromZerosAndPoles
  "  Transform a TransferFunction into a StateSpace description"
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Math.Complex;

  input ZerosAndPoles zp= Modelica_LinearSystems2.ZerosAndPoles({2+0*j}, {1+0*j,2+3*j,2-3*j}, 4);
protected
 input Complex j = Modelica_LinearSystems2.Math.Complex.j();
public
  output Boolean ok;

protected
  StateSpace ss=StateSpace(zp);

  ZerosAndPoles zp2;

algorithm
   Modelica.Utilities.Streams.print("zp = " + String(zp));
   Modelica.Utilities.Streams.print("ss = " + String(ss));
   ok := true;

end conversionFromZerosAndPoles;
