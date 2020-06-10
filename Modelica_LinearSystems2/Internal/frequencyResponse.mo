within Modelica_LinearSystems2.Internal;
encapsulated function frequencyResponse
  "Compute frequency response based on zeros and poles"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.Internal;
    import Modelica.Units.SI;

  input Real gain "Gain of transfer function";
  input Real Zeros[:,2]
    "Zeros as Real matrix (first column: real, second column imaginary values)";
  input Real Poles[:,2]
    "Poles as Real matrix (first column: real, second column imaginary values)";
  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
    "= true, if abszissa range is automatically determined";
  input SI.Frequency f_min(min=0) = 0.1
    "Minimum frequency value, if autoRange = false" annotation(Dialog(enable=not autoRange));
  input SI.Frequency f_max(min=0) = 10
    "Maximum frequency value, if autoRange = false" annotation(Dialog(enable=not autoRange));
  input Boolean Hz=true
    "= true, to compute abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
  input Boolean dB=false
    "= true, to compute magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));
  input Boolean logX=true
    "= true, to compute abszissa values for logarithmic scale"                       annotation(choices(checkBox=true));
  output Real f[nPoints] "Frequency vector (either in Hz or rad/s)";
  output Real A[nPoints] "Amplitude (either without unit or in dB)";
  output Modelica.Units.NonSI.Angle_deg phi[nPoints] "Angles in degree";
protected
  SI.AngularVelocity w[nPoints];
  SI.Angle phi_old;
  Integer info;
algorithm
  // Determine frequency vector f
  f := Internal.frequencyVector2(Zeros,Poles,nPoints,autoRange,
                                 f_min,f_max,logX);

  // Compute magnitude and phase at the frequency points
  phi_old := 0.0;
  for i in 1:nPoints loop
    w[i] := Modelica.Units.Conversions.from_Hz(f[i]);
    (A[i], phi_old, info) := Internal.frequencyEvaluate(gain, Zeros, Poles, 0, w[i]);
    phi[i] := Modelica.Units.Conversions.to_deg(phi_old);

    // Convert to other units, if required
    if not Hz then
       f[i] := w[i];
    end if;
    if dB then
       if A[i] <> 0 then
          A[i] := 20*log10(A[i]);
       else
          A[i] := -6000 "= 20*log10(1e-300)";
       end if;
    end if;
  end for;

  annotation (__Dymola_translate=true,Documentation(info="<html>
</html>"));
end frequencyResponse;
