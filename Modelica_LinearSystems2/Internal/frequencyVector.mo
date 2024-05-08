within Modelica_LinearSystems2.Internal;
function frequencyVector "Determine frequency vector for Bode plot"
  extends Modelica.Icons.Function;

  import Modelica.Units.SI;

  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
    "True, if abszissa range is automatically determined";
  input SI.Frequency f_min(min=0) = 0.1
    "Minimum frequency value, if autoRange = false"
    annotation (Dialog(enable=not autoRange));
  input SI.Frequency f_max(min=0) = 10
    "Maximum frequency value, if autoRange = false"
    annotation (Dialog(enable=not autoRange));
  input Complex numZeros[:]=fill(Complex(0), 0) "Zeros of numerator"
    annotation(Dialog(enable=autoRange));
  input Complex denZeros[:]=fill(Complex(0), 0) "Zeros of denominator"
    annotation(Dialog(enable=autoRange));
  input Boolean logX=true
    "=true: logarithmic scale; = false: linear scale of frequency vector";
  output SI.Frequency f[nPoints] "Frequency vector (automatic or manual)";
protected
  Real w_min;
  Real w_max;
  Real f_min2;
  Real f_min2b;
  Real f_max2;
  Real f_log[nPoints];
algorithm
  // Determine f_min2, f_max2 (auto or manual)
  if autoRange then
    (w_min,w_max) := Internal.frequencyRangeBode(numZeros, denZeros);
    f_min2 := Modelica.Units.Conversions.to_Hz(w_min);
    f_max2 := Modelica.Units.Conversions.to_Hz(w_max);
  else
    f_min2 := f_min;
    f_max2 := f_max;
  end if;

  // Compute vector of frequency points
  if logX then
    if f_min2 <= 0 then
      // Change minimum frequency point
      f_min2b :=f_min2;
      f_min2  :=f_max2/10^5;
      assert(false,"Minimum frequency changed from " + String(f_min2b) +
             " to " + String(f_min2) + " because log(0) is not defined",
             AssertionLevel.warning);
    end if;
    f_log := linspace(
      Modelica.Math.log10(f_min2),
      Modelica.Math.log10(f_max2),
      nPoints);
    for i in 1:nPoints loop
      f[i] := 10^f_log[i];
    end for;
  else
    f :=linspace(f_min,f_max,nPoints);
  end if;
end frequencyVector;
