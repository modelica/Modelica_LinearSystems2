within Modelica_LinearSystems2.Controller.Internal;
function FIR_coefficients "Calculates the FIR-filter coefficient vector"
  import Modelica_LinearSystems2.Utilities.Types.FilterType;
  import Modelica_LinearSystems2.Controller.Types.FIRspec;

  input FIRspec specType=FIRspec.MeanValue
    "Specification type of FIR filter";
  input Integer L(min=2) = 2 "Length of mean value filter" annotation(Dialog(enable=specType==FIRspec.MeanValue));
  input FilterType filterType=FilterType.LowPass "Type of filter" annotation (Dialog(enable=specType == FIRspec.Window));
  input Integer order(min=1) = 2 "Order of filter" annotation(Dialog(enable=specType==FIRspec.Window));
  input Modelica.Units.SI.Frequency f_cut=1 "Cut-off frequency"
    annotation (Dialog(enable=specType == FIRspec.Window));
  input Modelica.Units.SI.Time Ts(min=0) "Sampling time";
  input Types.Window window=Modelica_LinearSystems2.Controller.Types.Window.Rectangle
    "Type of window" annotation(Dialog(enable=specType==FIRspec.Window));
  input Real beta=2.12 "Beta-Parameter for Kaiser-window"
    annotation(Dialog(enable=specType==FIRspec.Window and window==Modelica_LinearSystems2.Controller.Types.Window.Kaiser));
  input Real a_desired[:]={1,1} "FIR filter coefficients" annotation(Dialog(enable=specType==FIRspec.Coefficients));
  output Real a[if specType==FIRspec.MeanValue then L else
                     (if specType == FIRspec.Window then if mod(order,2)>0 and filterType == FilterType.HighPass then order+2 else order+1 else
                     size(a_desired,1))] "Filter coefficients";

protected
  constant Real pi=Modelica.Constants.pi;
  Boolean isEven=mod(order,2)==0;
  Integer order2 = if not isEven and filterType == FilterType.HighPass then order+1 else order;
  Real Wc=2*pi*f_cut*Ts;
  Integer i;
  Real w[order2 + 1];
  Real k;
algorithm
 assert(f_cut<=1/(2*Ts),"The cut-off frequency f_cut may not be greater than half the sample frequency (Nyquist frequency), i.e. f_cut <= " + String(1/(2*Ts)) + " but is "+String(f_cut));
  if specType == FIRspec.MeanValue then
    a := fill(1/L, L);
  elseif specType == FIRspec.Window then
    w := Internal.FIR_window(order2 + 1, window, beta);
    for i in 1:order2 + 1 loop
      k := i - 1 - order2/2;
      if i - 1 == order2/2 then
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

  if not isEven and filterType == FilterType.HighPass then
    Modelica.Utilities.Streams.print("The requested order of the FIR filter in FIR_coefficients is odd and has been increased by one to get an even order filter\n");
  end if;

  annotation (
    Documentation(info="<html>
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
</p>
<pre>
     y(k) = a0*u(k) + a1*u(k-1) + a2*u(k-2) + ... + an*u(k-n)
</pre>
<p>
implies that the function outputs n+1 coefficients for a n-th order filter. The
coefficients can be weightened with different kind of windows: Rectangle, Bartlett,
Hann, Hamming, Blackman, Kaiser. The beta parameter is only needed by the Kaiser window.
</p>
</html>"));
end FIR_coefficients;
