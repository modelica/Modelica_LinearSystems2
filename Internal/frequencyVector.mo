within Modelica_LinearSystems2.Internal;
function frequencyVector "Determine frequency vector for Bode plot"
  import Modelica;
  import Modelica_LinearSystems2.Math.Complex;
  import LinearSystems = Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Internal;
  import SI = Modelica.SIunits;

  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
    "= true, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min(min=0) = 0.1
    "Minimum frequency value, if autoRange = false"                                                 annotation(Dialog(enable=not autoRange));
  input Modelica.SIunits.Frequency f_max(min=0) = 10
    "Maximum frequency value, if autoRange = false"                                                annotation(Dialog(enable=not autoRange));
  input Complex numZeros[:]=fill(Complex(0), 0) "Zeros of numerator" 
                                                                    annotation(Dialog(enable=autoRange));
  input Complex denZeros[:]=fill(Complex(0), 0) "Zeros of denominator" 
                                                                      annotation(Dialog(enable=autoRange));
  output SI.Frequency f[nPoints] "Frequency vector (automatic or manual)";
protected
  Real w_min;
  Real w_max;
  Real f_min2;
  Real f_max2;
  Real f_log[nPoints];
algorithm
  // Determine f_min2, f_max2 (auto or manual)
  if autoRange then
    (w_min,w_max) := Internal.frequencyRangeBode(numZeros, denZeros);
    f_min2 := Internal.to_Hz(w_min);
    f_max2 := Internal.to_Hz(w_max);

  else
    f_min2 := f_min;
    f_max2 := f_max;
  end if;
  w_min := Internal.from_Hz(f_min2);
  w_max := Internal.from_Hz(f_max2);

  // Compute logarithmic vector of frequency points
  f_log := linspace(
    Modelica.Math.log10(f_min2),
    Modelica.Math.log10(f_max2),
    nPoints);
  for i in 1:nPoints loop
    f[i] := 10^f_log[i];
  end for;
end frequencyVector;
