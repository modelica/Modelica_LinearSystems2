within Modelica_LinearSystems2.Internal;
record Eigenvalue
  "Record containing a eigen value or a pair of conjugated complex pair, respectively and additionally characteristics of the eigenvalue(s)"
  import Complex;
  extends Modelica.Icons.Record;

  Complex ev;
  Boolean imag = false;
  Boolean isStable = false;
  Boolean isControllable = false;
  Boolean isStabilizable = false;
  Boolean isObservable = false;
  Boolean isDetectable = false;
  Real frequency;
  Real damping;
  Real timeConstant;

encapsulated function constructor "Default constructor for eigenvalue"
  import Modelica;
  import Complex;
  import Modelica_LinearSystems2.Internal;
  import Modelica_LinearSystems2.Internal.Eigenvalue;

  input Complex ev_in=Complex(0);
  input Boolean isControllable=false;
  input Boolean isObservable = false;
  input Integer maxIndex1=0;
  input Integer maxIndex2=0;
  input Real Teps = 1e6
      "maximum time constant before regarded as infinity, i.e. a real eigenvalue is zero";

  output Eigenvalue ev;

  protected
  Boolean isStable = ev_in.re<0;
  Boolean isImag = abs(ev_in.im)>Modelica.Constants.eps;
  Real abs_ev = (ev_in.re^2+ev_in.im^2)^0.5;

algorithm
  ev.ev := ev_in;
  ev.imag := isImag;
  ev.isStable := isStable;
  ev.isControllable := isControllable;
  ev.isStabilizable := isStable or isControllable;
  ev.isObservable := isObservable;
  ev.isDetectable := isStable or isObservable;
  ev.frequency := if ev.imag then abs_ev/(2*Modelica.Constants.pi) else 0.0;
  ev.damping := if ev.imag then if abs_ev>Modelica.Constants.eps then -ev.ev.re/abs_ev else 0.0 else 1.0;
  ev.timeConstant := if ev.imag then 0.0 else if abs(ev.ev.re) > 1/Teps then 1/abs(ev.ev.re) else Teps;
end constructor;

end Eigenvalue;
