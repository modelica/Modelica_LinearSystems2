within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode1
  "Example for construction of a ZerosAndPoles system and plot of the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input Modelica.SIunits.Frequency f_cut=100 "PT1 with cut-off frequency f_cut";
  output Boolean ok;
protected
  Modelica.SIunits.AngularVelocity w=2*Modelica.Constants.pi*f_cut;
  TransferFunction tf=TransferFunction(n={w}, d={1,w});
  ZerosAndPoles zp=ZerosAndPoles(tf);
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(zp);

  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example shows how to construct a zeros and poles system and to plot the Bode diagram
with automatic determination of the frequency range to plot.
</p>
</html>
"));
end plotBode1;
