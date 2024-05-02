within Modelica_LinearSystems2.Internal;
function frequencyRangeBode2
  "Determine min. and max. frequencies for Bode plot (Zeros and Poles as Real matrix)"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Internal;

  input Real Zeros[:,2]
    "Zeros as Real matrix (first column: real, second column imaginary values)";
  input Real Poles[:,2]
    "Poles as Real matrix (first column: real, second column imaginary values)";
  output Modelica.Units.SI.AngularVelocity w_min "Minimum angular frequency";
  output Modelica.Units.SI.AngularVelocity w_max "Maximum angular frequency";
protected
  Real phi_min=Modelica.Units.Conversions.from_deg(5);
  Real real_min=1.0e-4;
  Real pi=Modelica.Constants.pi;
  Integer n_num;
  Integer n_den;
  Real w_min1;
  Real w_min2;
  Real w_max1;
  Real w_max2;
  Boolean useFullRange1;
  Boolean useFullRange2;
algorithm
  // Compute frequencies for numerator
  n_num := size(Zeros, 1);
  if n_num > 0 then
    (w_min1,w_max1,useFullRange1) := Internal.frequencyRangeZeros2(
                                        Zeros, phi_min, real_min);
  end if;

  // Compute frequencies for denominator
  n_den := size(Poles, 1);
  if n_den > 0 then
    (w_min2,w_max2,useFullRange2) := Internal.frequencyRangeZeros2(
                                       Poles, phi_min, real_min);
  end if;

  // Use largest range
  if n_num == 0 and n_den == 0 then
    w_min :=Modelica.Units.Conversions.from_Hz(0.1);
    w_max :=Modelica.Units.Conversions.from_Hz(1);
  elseif n_num == 0 then
    w_min := w_min2;
    w_max := w_max2;
  elseif n_den == 0 then
    w_min := w_min1;
    w_max := w_max1;
  else
    if useFullRange1 and useFullRange2 then
      w_min := min(w_min1, w_min2);
      w_max := max(w_max1, w_max2);
    elseif useFullRange1 then
      w_min :=w_min1;
      w_max :=w_max1;
    else
      w_min :=w_min2;
      w_max :=w_max2;
    end if;
  end if;
end frequencyRangeBode2;
