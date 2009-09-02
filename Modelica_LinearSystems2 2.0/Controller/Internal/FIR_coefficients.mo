within Modelica_LinearSystems2.Controller.Internal;
function FIR_coefficients "Calculates the FIR-filter coefficient vector"
  import Modelica_LinearSystems2.Types.FilterType;
  import Modelica_LinearSystems2.Controller.Types.FIRspec;

  input FIRspec specType=Modelica_LinearSystems2.Controller.Types.FIRspec.MeanValue
    "Specification type of FIR filter";
  input Integer L(min=2) = 2 "Length of mean value filter" annotation(Dialog(enable=specType==FIRspec.MeanValue));
  input Modelica_LinearSystems2.Types.FilterType filterType=
      Modelica_LinearSystems2.Types.FilterType.LowPass "Type of filter"
                            annotation(Dialog(enable=specType==FIRspec.Window));
  input Integer order(min=1) = 2 "Order of filter" annotation(Dialog(enable=specType==FIRspec.Window));
  input Modelica.SIunits.Frequency f_cut=1 "Cut-off frequency" annotation(Dialog(enable=specType==FIRspec.Window));
  input Modelica.SIunits.Time Ts(min=0) "Sampling time";
  input Types.Window window=Modelica_LinearSystems2.Controller.Types.Window.Rectangle
    "Type of window" annotation(Dialog(enable=specType==FIRspec.Window));
  input Real beta=2.12 "Beta-Parameter for Kaiser-window"
    annotation(Dialog(enable=specType==FIRspec.Window and window==Modelica_LinearSystems2.Controller.Types.Window.Kaiser));
  input Real a_desired[:]={1,1} "FIR filter coefficients" annotation(Dialog(enable=specType==FIRspec.Coefficients));
  output Real a[if specType==FIRspec.MeanValue then L else
                     (if specType == FIRspec.Window then order+1 else
                     size(a_desired,1))] "Filter coefficients";

  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Window(
      x=0.22,
      y=0.24,
      width=0.64,
      height=0.61),
    Documentation(info="<HTML>
<p>
The FIR-filter synthesis based on the window method. The coefficients are
calculated through a fourier series approximation of the desired amplitude
characteristic. Due to the fact that the Fourier series is truncated, there
will be discontinuities in the magnitude of the filter. Especial at the edge
of the filter the ripple is concentrated (Gibbs-effect). To counteract this,
the filter coefficients are convolved in the frequency domain with the spectrum
of a window function, thus smoothing the edge transitions at any discontinuity.
This convolution in the frequency domain is equivalent to multiplying the filter
coefficients with the window coefficients in the time domain.
</p>
<p>
The filter equation
<pre>
     y(k) = a0*u(k) + a1*u(k-1) + a2*u(k-2) + ... + an*u(k-n)
</pre>
implies that the function outputs n+1 coefficients for a n-th order filter. The
coefficients can be weightened with different kind of windows: Rectangle, Bartlett,
Hann, Hamming, Blackman, Kaiser <br>
The beta parameter is only needed by the Kaiser window.
</p>
</HTML>
"));
protected
  constant Real pi=Modelica.Constants.pi;
  Real Wc=2*pi*f_cut*Ts;
  Integer i;
  Real w[order + 1];
  Real k;
algorithm
  if specType == FIRspec.MeanValue then
     a := fill(1/L, L);
  elseif specType == FIRspec.Window then
     Modelica.Utilities.Streams.error(
         "  There seems to be a bug when calculating the FIR coefficients for\n"+
         "  specType = FIRspec.Window. Therefore, the calculation is temporarily\n" +
         "  switched off");
     (w) := Internal.FIR_window(order + 1, window, beta);
     for i in 1:order + 1 loop
       k := i - 1 - order/2;
       if i - 1 == order/2 then
         a[i] := if filterType == FilterType.LowPass then Wc*w[i]/pi else
                 w[i] - Wc*w[i]/pi;
       else
         a[i] := if filterType == FilterType.LowPass then sin(k*Wc)*
           w[i]/(k*pi) else w[i]*(sin(k*pi) - sin(k*Wc))/(k*pi);
       end if;
     end for;
  else
     a := a_desired;
  end if;
end FIR_coefficients;
