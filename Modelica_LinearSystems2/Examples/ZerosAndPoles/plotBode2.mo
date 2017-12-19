within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode2 "Example for a Bode plot of PT2 transfer function with zero damping"
  extends Modelica.Icons.Function;

  import Complex;
  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input Modelica.SIunits.Frequency f_cut=10
    "Cut-off frequency of denominator PT2";
  input Real D=0 "Damping of denominator PT2";
  input Real k=1 "Gain";
  input Integer nPoints=1000;
  output Boolean ok;
protected
  Modelica.SIunits.AngularVelocity w=2*Modelica.Constants.pi*f_cut;
  TransferFunction tf=TransferFunction(n={k*w^2}, d={1,2*D*w,w*w});
  ZerosAndPoles zp=ZerosAndPoles(tf);
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(zp, nPoints);

  ok := true;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example shows how to construct a PT2 system with zero damping and to plot the Bode diagram.
</p>
</html>
"));
end plotBode2;
