within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function conv_zp2ss

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input ZerosAndPoles zp;
  input Integer n=0;

protected
  StateSpace ss1=StateSpace(zp);
  TransferFunction tf1=StateSpace.Conversion.toTransferFunction(ss1);
  TransferFunction tf2=ZerosAndPoles.Conversion.toTransferFunction(zp);
  ZerosAndPoles zp2=ZerosAndPoles(tf2);
  StateSpace ss2=StateSpace(zp2);
  TransferFunction tf3=StateSpace.Conversion.toTransferFunction(ss2);

algorithm
  if n>0 then
  Modelica.Utilities.Streams.print("\nNr. " + String(n));
  end if;

  Modelica.Utilities.Streams.print("\nzp = " + String(zp));
  Modelica.Utilities.Streams.print("tf1 = " + String(tf1));
  Modelica.Utilities.Streams.print("tf2 = " + String(tf2));
  Modelica.Utilities.Streams.print("\nzp2 = " + String(zp2));
  Modelica.Utilities.Streams.print("tf3 = " + String(tf3));
  Modelica.Utilities.Streams.print("tf2 = " + String(tf2));

end conv_zp2ss;
