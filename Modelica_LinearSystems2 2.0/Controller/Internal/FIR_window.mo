within Modelica_LinearSystems2.Controller.Internal;
function FIR_window "Calculation of n-point weighting window for FIR filter"

  import Modelica_LinearSystems2.Controller.Types.Window;
  input Integer L "Number of Points";
  input Window window "Type of window";
  input Real beta=2.12 "Beta-Parameter for Kaiser-window";
  output Real a[L] "output vector";
protected
  Integer i=0;
  constant Real pi=Modelica.Constants.pi;
  Real k;
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Window(
      x=0.28,
      y=0.15,
      width=0.64,
      height=0.71),
    Documentation(info="<HTML>
<p>
Weighting windows are used for digital filter design or spectrum estimation (e.g. DFT)
to increase the quality. In designing FIR-Filter the main role of windowing is to remove
non-ideal effects caused by the endless number of filter coefficients (Gibbs phenomenon).
Multiplying the coefficients with a window damps the coefficients at the beginning and at
the end.
</p>
<p>
The function outputs a L-point vector for a given kind of window. The parameter \"beta\" is
only needed by the Kaiser window. The types of windows are:
</p>
<OL>
<LI>Rectangle</LI>
<LI>Bartlett</LI>
<LI>Hann</LI>
<LI>Hamming</LI>
<LI>Blackman</LI>
<LI>Kaiser</LI>
</OL>
</HTML>
"));
algorithm
  if window <> Window.Rectangle then
    for i in 1:L loop
      k := i - 1 - (L - 1)/2;
      if window == Window.Bartlett then
        a[i] := 1 - 2*abs(k)/(L - 1);
      elseif window == Window.Hann then
        a[i] := 0.5 + 0.5*cos(2*pi*k/(L - 1));
      elseif window == Window.Hamming then
        a[i] := 0.54 + 0.46*cos(2*pi*k/(L - 1));
      elseif window == Window.Blackman then
        a[i] := 0.42 + 0.5*cos(2*pi*k/(L - 1)) + 0.08*cos(4*pi*k/(L - 1));
      elseif window == Window.Kaiser then
        k := 2*beta*sqrt((i - 1)*(L - i))/(L - 1);
        a[i] := Internal.bessel0(k)/ Internal.bessel0(beta);
      else
        Modelica.Utilities.Streams.error("window = " + String(window) + " not known");
      end if;
    end for;
  else
    a := ones(L);
  end if;
end FIR_window;
