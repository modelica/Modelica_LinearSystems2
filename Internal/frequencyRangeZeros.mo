within Modelica_LinearSystems2.Internal;
function frequencyRangeZeros
  "Determine min. and max. frequencies for a vector of zeros (numerator or denominator zeros)"
  import Modelica;
  import Modelica.Math;
  import SI = Modelica.SIunits;
  import Modelica_LinearSystems2.Math.Complex;

  input Complex z[:] "Vector of zeros";
  input SI.Angle phi_min(min=10*Modelica.Constants.eps)=
    Modelica.SIunits.Conversions.from_deg(5) "Minimum phase angle";
  input Real real_min(min=0) = 1.e-4 "[r| < real_min are treated as |real_min|";
  output Modelica.SIunits.AngularVelocity w_min "Minimum frequency";
  output Modelica.SIunits.AngularVelocity w_max "Maximum frequency";

protected
  Integer nz=size(z, 1);
  Real tan_min;
  Real tan_max1;
  Real tan_max2;
  Real z_re;
  Real w_min1;
  Real w_max1;
  Real z_abs2;
  Real k1;
  Real k2;
algorithm
  /* - Real zero:
       tan(phi_desired) = w/re -> w = re*tan(phi_desired)
  
     - Conjugate complex zero:
       tan(phi1) = (w + im)/re
       tan(phi2) = (w - im)/re
       phi1 + phi2 = phi_desired -> tan(phi1 + phi2) = (tan(phi1) + tan(phi2))/(1-tan(phi1)*tan(phi2))
       -> after a longer derivation it follows:
          w = sqrt( re^2/tan(phi_desired)^2 + (re^2 + im^2) ) - re/tan(phi_desired)
  */
  assert(nz > 0, "Vector z of zeros has dimension 0, This is not allowed");
  tan_min := Math.tan(phi_min);
  tan_max1 := Math.tan(Modelica.Constants.pi/2 - phi_min);
  tan_max2 := Math.tan(Modelica.Constants.pi - phi_min);
  for i in 1:size(z, 1) loop
    if z[i].im >= 0.0 then
      z_re := max(abs(z[i].re), real_min);
      if z[i].im > 0.0 then
        z_abs2 := z_re^2 + z[i].im^2;
        k1 := z_re/tan_min;
        w_min1 := sqrt(k1^2 + z_abs2) - k1;
        k2 := z_re/tan_max2;
        w_max1 := sqrt(k2^2 + z_abs2) - k2;
      else
        w_min1 := z_re*tan_min;
        w_max1 := z_re*tan_max1;
      end if;

      if i == 1 then
        w_min := w_min1;
        w_max := w_max1;
      else
        w_min := min(w_min, w_min1);
        w_max := max(w_max, w_max1);
      end if;
    end if;
  end for;
  annotation (Documentation(info="<html>
<p>
This function estimates a useful frequency range for the
Bode plot of a vector of zeros (numerator or denominator zeros).
This frequency range is estimated such that the phase angle 
of <b>one</b> zero is in the range:
<p>
<pre>
   phi_min/n_zeros &le; |phase angle| &le; pi/2 - phi_min/n_zeros
</pre>
<p>
where n_zeros is the number of zeros.
Note, the phase angle of one zero for a frequency of 0 up to infinity
is in the range:
</p>
<pre>
   0 &le; |phase angle| &le; pi/2
</pre>
<p>
Therefore, the frequency range is estimated
such that the essential part of the phase angle (defined by phi_min)
is present.
</p>
<p>
If the real part of a complex zero vanishes 
(i.e., the zero is located on the imaginary axis), 
the maximum value of the bode plot magnitude of the zero 
is infinity. In order to avoid difficulties, zeros close to
the imaginary axis are shifted by the input argument
real_min along the real axis.
</p>
</html>"));
end frequencyRangeZeros;
