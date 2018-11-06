within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode3
  "Example for construction of a ZerosAndPoles system and plot of the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Polynomial;

  input Modelica.SIunits.Frequency f_cut_num=100
    "Cut-off frequency of numerator PT2 and PT1";
  input Real D_num=0.1 "Damping of numerator PT2";
  input Modelica.SIunits.Frequency f_cut_den=10
    "Cut-off frequency of denominator PT2 and PT1";
  input Real D_den=0.1 "Damping of denominator PT2";
  input Real k=1 "Gain";
  output Boolean ok;
protected
  Modelica.SIunits.AngularVelocity w1=2*Modelica.Constants.pi*f_cut_num;
  Modelica.SIunits.AngularVelocity w2=2*Modelica.Constants.pi*f_cut_den;

  Polynomial pn1=Polynomial(k/w1^3*{1,w1});
  Polynomial pn2=Polynomial({1,2*D_num*w1,w1*w1});
  Polynomial pn3=pn1*pn2;
  Polynomial pd1=Polynomial({1,w2});
  Polynomial pd2=Polynomial({1,2*D_den*w2,w2*w2});
  Polynomial pd3=pd1*pd2;

  TransferFunction tf1=TransferFunction(n=pn3, d=Polynomial({1}));
  TransferFunction tf2=TransferFunction(n=Polynomial({k*w2^3}), d=pd3);
  TransferFunction tf3=tf1*tf2;

  ZerosAndPoles zp1=ZerosAndPoles(tf1);
  ZerosAndPoles zp2=ZerosAndPoles(tf2);
  ZerosAndPoles zp3=ZerosAndPoles(tf3);

  Modelica.SIunits.Frequency f_min=2;
  Modelica.SIunits.Frequency f_max=200;
algorithm
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    zp1,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    zp2,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(
    zp3,
    autoRange=false,
    f_min=f_min,
    f_max=f_max);

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows how to construct a zeros and poles system and to plot the Bode diagram
with automatic determination of the frequency range to plot.
</p>
</html>
"));
end plotBode3;
