within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotPolesAndZeros
  "Example for plotting poles and zeros of a ZerosAndPoles transfer function"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

protected
  TransferFunction s = TransferFunction.s();
  TransferFunction tf = (s^3 + 4*s + 1)/(s^4 + 2*s^3 + 3*s^2 + 4*s);
  ZerosAndPoles zp = ZerosAndPoles(tf);
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.polesAndZeros(
    zp=zp,
    defaultDiagram = Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
      heading="Poles and zeros of " + String(tf)),
      device=Modelica_LinearSystems2.Utilities.Plot.Records.Device(xTopLeft=50, yTopLeft=30));

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This example shows how to plot a pole-zero-map of a transfer function in ZerosAndPoles representation given internally.
</p>
</html>
"));
end plotPolesAndZeros;
