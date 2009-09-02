within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode1
  "Construct a ZerosAndPoles system and plot the Bode diagram with automatic determination of the frequency range to plot"
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica.SIunits.Frequency f_cut=100 "PT1 with cut-off frequency f_cut";
  output Boolean ok;
  annotation (interactive=true);
protected
  Modelica.SIunits.AngularVelocity w=2*Modelica.Constants.pi*f_cut;
  TransferFunction tf=TransferFunction(n={w}, d={1,w});
  ZerosAndPoles zp=ZerosAndPoles(tf);
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(    zp);
  ok := true;

equation

end plotBode1;
