within Modelica_LinearSystems2.Internal;
function frequencyRangeZeros2
  "Determine min. and max. frequencies for a vector of zeros (numerator or denominator zeros; provided as Real matrix)"
  import Modelica;
  import Modelica.Math;
  import      Modelica.Units.SI;

  input Real Zeros[:,2]
    "Zeros as Real matrix (first column: real, second column imaginary values)";
  input SI.Angle phi_min(min=10*Modelica.Constants.eps)=
    Modelica.Units.Conversions.from_deg(5) "Minimum phase angle";
  input Real real_min(min=0) = 1.e-4 "|r| < real_min are treated as |real_min|";
  output Modelica.Units.SI.AngularVelocity w_min "Minimum frequency";
  output Modelica.Units.SI.AngularVelocity w_max "Maximum frequency";
  output Boolean useFullRange = true;
protected
  Integer nz=size(Zeros, 1);
  Real tan_min1;
  Real tan_min2;
  Real tan_max1;
  Real tan_max2;
  Real z_re;
  Real w_min1;
  Real w_max1;
  Real z_abs2;
  Real k1;
  Real k2;
  Boolean first = true;
algorithm
  /* - Real zero:
       tan(phi_desired) = w/|re| -> w = |re|*tan(phi_desired)

     - Conjugate complex zero:
       tan(phi1) = (w + im)/|re|
       tan(phi2) = (w - im)/|re|
       phi1 + phi2 = phi_desired ->
          tan(phi1 + phi2) = (tan(phi1) + tan(phi2))/(1-tan(phi1)*tan(phi2))
                           = ((w + im)/re + (w - im)/re) / ( 1 - (w+im)/re*(w-im)/re )
                           = 2*w/re / (re^2 - (w^2 - im^2))/re^2
                           = 2*w*re / (re^2 - w^2 + im^2)
          tan(phi_desired)*(re^2 + im^2 - w^2) = 2*w*re
          tan(phi_desired)*w^2 + 2*w*re - tan(phi_desired)*(re^2 + im^2) = 0
          w^2 + (2*re/tan(phi_desired))*w - (re^2 + im^2) = 0
          w = -re/tan(phi_desired) + sqrt(re^2/tan(phi_desired)^2 + re^2 + im^2)
  */
  assert(nz > 0, "Matrix Zeros has dimension 0, This is not allowed");
  tan_min1 := Math.tan(phi_min);
  tan_min2 := Math.tan(2*phi_min);
  tan_max1 := Math.tan(Modelica.Constants.pi/2 - phi_min);
  tan_max2 := Math.tan(Modelica.Constants.pi - 2*phi_min);
  for i in 1:nz loop
    if Zeros[i,2] >= 0.0 and abs(Zeros[i,1]) >= real_min then
      //z_re := if abs(z[i].re) <= real_min then -real_min else z[i].re;
      z_re :=max(abs(Zeros[i,1]), real_min);
      if Zeros[i,2] > 0.0 then
        z_abs2 := z_re^2 + Zeros[i,2]^2;
        k1 := z_re/tan_min2;
        w_min1 := sqrt(k1^2 + z_abs2) - k1;
        k2 := z_re/tan_max2;
        w_max1 := sqrt(k2^2 + z_abs2) - k2;
        /*
        phi1 :=Modelica.SIunits.Conversions.to_deg(Modelica.Math.atan((w_max1 + z[i].im)/z_re));
        phi2 :=Modelica.SIunits.Conversions.to_deg(Modelica.Math.atan((w_max1 - z[i].im)/z_re));
        Modelica.Utilities.Streams.print("phi1, phi2, 180-sum = " + String(phi1) + ", " + String(phi2)+", " +String(180-phi1-phi2));
        */
      else
        w_min1 := z_re*tan_min1;
        w_max1 := z_re*tan_max1;
      end if;

      if first then
        first :=false;
        w_min := w_min1;
        w_max := w_max1;
      else
        w_min := min(w_min, w_min1);
        w_max := max(w_max, w_max1);
      end if;
    end if;
  end for;

  if first then
     useFullRange := false;
     w_min :=Modelica.Units.Conversions.from_Hz(0.1);
     w_max :=Modelica.Units.Conversions.from_Hz(1);
  end if;

  annotation (Documentation(info="<html>
<p>
This function estimates a useful frequency range for the
Bode plot of a vector of zeros (numerator or denominator zeros).
This frequency range is estimated such that the phase angle
of <b>one</b> zero is in the range:
</p>
<blockquote><pre>
phi_min/n_zeros &le; |phase angle| &le; pi/2 - phi_min/n_zeros
</pre></blockquote>

<p>
where n_zeros is the number of zeros.
Note, the phase angle of one zero for a frequency of 0 up to infinity
is in the range:
</p>
<blockquote><pre>
0 &le; |phase angle| &le; pi/2
</pre></blockquote>

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
end frequencyRangeZeros2;
