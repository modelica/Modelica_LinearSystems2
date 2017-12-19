within Modelica_LinearSystems2.Examples.TransferFunction;
function plotPolesAndZeros
  "Example for plotting poles and zeros of two transfer functions"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica.ComplexMath.j;

protected
  TransferFunction s = TransferFunction.s();
  TransferFunction tf1 = (s-2)*(s+4)/( (s+1)*(s+3)*TransferFunction({-0.5+2*j, -0.5-2*j}));
  TransferFunction tf2 = (s^3 + 4*s + 1)/(s^4 + 2*s^3 + 3*s^2 + 4*s);

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.polesAndZeros(
    tf=tf1,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
      heading="Poles and zeros of " + String(tf1)));

  /*
  Modelica_LinearSystems2.TransferFunction.Plot.polesAndZeros(
    tf=tf2,
    defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
      heading="Poles and zeros of " + String(tf2)),
    device=Modelica_LinearSystems2.Utilities.Plot.Records.Device(xTopLeft=50, yTopLeft=30));
  */
  annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This example shows how to plot a pole-zero-map of a transfer function given internally.
</p>
</html>
"));
end plotPolesAndZeros;
