within Modelica_LinearSystems2.Internal;
function frequencyRangeBode "Determine min. and max. frequencies for Bode plot"
  import Modelica;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Internal;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica.Utilities.Streams;

  input Complex numZeros[:] "Zeros of numerator";
  input Complex denZeros[:] "Zeros of denominator";
  output Modelica.SIunits.AngularVelocity w_min;
  output Modelica.SIunits.AngularVelocity w_max;
protected
  Real phi_min=Modelica.SIunits.Conversions.from_deg(3);
  Real real_min=1.0e-4;
  Real pi=Modelica.Constants.pi;
  Integer n_num;
  Integer n_den;
  Real w_min1;
  Real w_min2;
  Real w_max1;
  Real w_max2;
  Real f_min;
  Real f_max;
algorithm
  // Compute frequencies for numerator
  n_num := size(numZeros, 1);
  if n_num > 0 then
    (w_min1,w_max1) := Internal.frequencyRangeZeros(
      numZeros,
      phi_min,
      real_min);
  end if;

  // Compute frequencies for denominator
  n_den := size(denZeros, 1);
  if n_den > 0 then
    (w_min2,w_max2) := Internal.frequencyRangeZeros(
      denZeros,
      phi_min,
      real_min);
  end if;

  // Use largest range
  if n_num == 0 and n_den == 0 then
    w_min := 0.1;
    w_max := 10;
  elseif n_num == 0 then
    w_min := w_min2;
    w_max := w_max2;
  elseif n_den == 0 then
    w_min := w_min1;
    w_max := w_max1;
  else
    w_min := min(w_min1, w_min2);
    w_max := max(w_max1, w_max2);
  end if;
  f_min := to_Hz(w_min);
  f_max := to_Hz(w_max);
end frequencyRangeBode;
