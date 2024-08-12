within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotBode3
  "Construct three zeros and poles systems and plot the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.Math.Polynomial;

  input Modelica.Units.SI.Frequency f_cut_num=100 "Cut-off frequency of numerator PT2 and PT1";
  input Real D_num=0.1 "Damping of numerator PT2";
  input Modelica.Units.SI.Frequency f_cut_den=10 "Cut-off frequency of denominator PT2 and PT1";
  input Real D_den=0.1 "Damping of denominator PT2";
  input Real k=1 "Gain";
  output Boolean ok;
protected
  Modelica.Units.SI.AngularVelocity w1=2*Modelica.Constants.pi*f_cut_num;
  Modelica.Units.SI.AngularVelocity w2=2*Modelica.Constants.pi*f_cut_den;

  Polynomial p1=Polynomial(k/w1^3*{1,w1});
  Polynomial p2=Polynomial({1,2*D_num*w1,w1*w1});
  Polynomial pn1=p1*p2;
  Polynomial p3=Polynomial({1,w2});
  Polynomial p4=Polynomial({1,2*D_den*w2,w2*w2});
  Polynomial pd2=p3*p4;

  TransferFunction tf1=TransferFunction(n=pn1, d=Polynomial({1}));
  TransferFunction tf2=TransferFunction(n=Polynomial({k*w2^3}), d=pd2);
  TransferFunction tf3=tf1*tf2;

  ZerosAndPoles zp1=ZerosAndPoles(tf1);
  ZerosAndPoles zp2=ZerosAndPoles(tf2);
  ZerosAndPoles zp3=ZerosAndPoles(tf3);

  Modelica.Units.SI.Frequency f_min=2;
  Modelica.Units.SI.Frequency f_max=200;
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

  annotation (__Dymola_interactive=true);
end plotBode3;
