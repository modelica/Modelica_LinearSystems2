within Modelica_LinearSystems2.Internal;
record PolesAndZeros
  "Record containing poles and zeros of a system in two real vectors containing the real parts and the imaginary parts respctively"
  import Modelica_LinearSystems2.Math.Complex;
  extends Modelica.Icons.Record;

  Real p_real[:];
  Real p_im[:];
  Real z_real[:];
  Real z_im[:];

  Integer norz_p "number of real zeros of p";
  Integer norz_z "number of real zeros of z";
  encapsulated function constructor "Default constructor for poles and zeros"
    import Modelica;
    import Modelica_LinearSystems2.Math.Complex;
    import PolesAndZeros2 = Modelica_LinearSystems2.Internal.PolesAndZeros;
    import Modelica_LinearSystems2.Internal;

    input Complex p[:]={Complex(0, 0)} "Complex zeros";
    input Complex z[:]={Complex(0, 0)} "Complex zeros";
    output PolesAndZeros2 pz(
      redeclare Real p_real[size(p, 1)],
      redeclare Real p_im[size(p, 1)],
      redeclare Real z_real[size(z, 1)],
      redeclare Real z_im[size(z, 1)]);
  algorithm
    pz.p_real := p[:].re;
    pz.p_im := p[:].im;
    pz.z_real := z[:].re;
    pz.z_im := z[:].im;
    pz.norz_p := Internal.numberOfRealZeros(p);
    pz.norz_z := Internal.numberOfRealZeros(z);

  end constructor;

end PolesAndZeros;
