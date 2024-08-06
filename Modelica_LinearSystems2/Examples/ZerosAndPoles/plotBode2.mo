within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode2 "Bode plot of PT2 transfer function with zero damping"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input Modelica.Units.SI.Frequency f_cut=10
    "Cut-off frequency of denominator PT2";
  input Real D=0 "Damping of denominator PT2";
  input Real k=1 "Gain";
  input Integer nPoints=1000;
  output Boolean ok;
protected
  Modelica.Units.SI.AngularVelocity w=2*Modelica.Constants.pi*f_cut;
  TransferFunction tf=TransferFunction(n={k*w^2}, d={1,2*D*w,w*w});
  ZerosAndPoles zp=ZerosAndPoles(tf);
algorithm
  ZerosAndPoles.Plot.bode(zp, nPoints);

  ok := true;

  annotation (__Dymola_interactive=true);
end plotBode2;
