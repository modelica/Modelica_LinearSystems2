within Modelica_LinearSystems2.WorkInProgress;
record DiscreteZerosAndPoles
  "Discrete zeros and poles description of a single input, single output system (data + operations)"
  extends Modelica.Icons.Record;

  Real k=1.0 "Multiplicative factor of transfer function"
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real n1[:] "[p^0] coefficients of 1st order numerator polynomials"
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real n2[:,2] "[p,p^0] coefficients of 2nd order numerator polynomials"
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real d1[:] "[p^0] coefficients of 1st order denominator polynomials"
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real d2[:,2] "[p,p^0] coefficients of 2nd order denominator polynomials"
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));

    Modelica.SIunits.Time Ts "Sample time"
       annotation(Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

/* If the numerator polynomial has no coefficients, the transfer function
   is zero. The denominator polynomial must always have at
   least one coefficient, such as {1}
*/

  String uName="u" "Name of input signal"    annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal"  annotation(Dialog(group="Signal names"));

  encapsulated operator 'constructor' "Generate a ZerosAndPoles object"
    import Modelica_LinearSystems2;

    encapsulated function fromReal
      "Generate a ZerosAndPoles transfer function from a Real value"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteZerosAndPoles dzp(
        redeclare Real n1[0],
        redeclare Real n2[0,2],
        redeclare Real d1[0],
        redeclare Real d2[0,2]) "= r";

    algorithm
      dzp.k := r;
      dzp.Ts := Ts;
      dzp.method := method;
      dzp.uName := uName;
      dzp.yName := yName;

    end fromReal;

  encapsulated function fromZerosAndPoles
      "Generate a ZerosAndPoles transfer function from a set of zeros and poles"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Internal;
    import Modelica.Utilities.Streams.print;

    input Complex z[:] = fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
    input Complex p[:] = fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
    input Real k = 1.0 "Constant multiplied with transfer function";
    input Modelica.SIunits.Time Ts = 1 "Sample time";
    input Modelica_LinearSystems2.Utilities.Types.Method method=
      Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
    input String uName = "" "Input name";
    input String yName = "" "Output name";
    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[Internal.numberOfRealZeros(z)],
      redeclare Real n2[integer((size(z, 1) - Internal.numberOfRealZeros(z))/2), 2],
      redeclare Real d1[Internal.numberOfRealZeros(p)],
      redeclare Real d2[integer((size(p, 1) - Internal.numberOfRealZeros(p))/2), 2])
      "ZerosAndPoles transfer functions of the zeros, poles and k";

    protected
    Integer n_n1=size(dzp.n1, 1);
    Integer n_d1=size(dzp.d1, 1);
    Integer n_n2=size(dzp.n2, 1);
    Integer n_d2=size(dzp.d2, 1);
    Integer i;
    Integer j;
    Complex z_reordered[size(z, 1)] "Reordered zeros";
    Complex p_reordered[size(p, 1)] "Reordered poles";
    Integer nz_real "Number of real zeros";
    Integer np_real "Number of real poles";
  algorithm
    //Extra input (string) added
    (z_reordered,nz_real) := Internal.reorderZeros(z, "");
    (p_reordered,np_real) := Internal.reorderZeros(p, "");

    // Numerator
    for i in 1:n_n1 loop
      dzp.n1[i] := -z_reordered[i].re;
    end for;

    j := 1;
    for i in n_n1 + 1:2:size(z, 1) loop
      dzp.n2[j, :] := {-2*z_reordered[i].re,z_reordered[i].re^2 + z_reordered[i].im
        ^2};
      j := j + 1;
    end for;

    // Denominator
    for i in 1:n_d1 loop
      dzp.d1[i] := -p_reordered[i].re;
    end for;

    j := 1;
    for i in n_d1 + 1:2:size(p, 1) loop
      dzp.d2[j, :] := {-2*p_reordered[i].re,p_reordered[i].re^2 + p_reordered[i].im
        ^2};
      j := j + 1;
    end for;

    dzp.Ts := Ts;
    dzp.method := method;
    dzp.k := k;
    dzp.uName := uName;
    dzp.yName := yName;

    annotation (Documentation(info="<html>
<p>
This function constructs a transfer function from denominator
and numerator zeros, as well as a gain.
Example:
</p>
<pre>                          (s+1)
  zp = 4* -------------------------------------
           (s - 1)*(s - (2+j*3))*(s - (2-j*3))
</pre>
<p>
with j=sqrt(-1), is defined as
</p>
<pre>
   <b>import</b> Modelica_LinearSystems2.Math.Complex;
   <b>import</b> Modelica_LinearSystems2.ZerosAndPoles;

   zp = ZerosAndPoles(z = {Complex(-1,0)},
                      p = {Complex(1,0),
                           Complex(2,3),
                           Complex(2,-3)},
                           k=4);
</pre>
<p>
Since only transfer functions with real coefficients are supported,
complex roots must be defined as conjugate complex pairs.
It is required that complex conjugate pairs must directly
follow each other as above. An error occurs if this is not the case.
</p>
</html>"));
  end fromZerosAndPoles;

    function fromDiscreteTransferFunction =
        Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction.Conversion.toDiscreteZerosAndPoles;
    encapsulated function fromFactorization
      "Generate a ZerosAndPoles object from first and second order polynomials"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

      input Real n1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order numerator polynomials"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real n2[:,2]=fill(
              0,
              0,
              2) "[p,p^0] coefficients of 2nd order numerator polynomials"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order denominator polynomials"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d2[:,2]=fill(
              0,
              0,
              2) "[p,p^0] coefficients of 2nd order denominator polynomials"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real k=1.0 "Multiplicative factor of transfer function"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteZerosAndPoles dzp(
        redeclare Real n1[size(n1, 1)],
        redeclare Real n2[size(n2, 1),2],
        redeclare Real d1[size(d1, 1)],
        redeclare Real d2[size(d2, 1),2]) "ZerosAndPoles transfer function";
    algorithm
      dzp.n1 := n1;
      dzp.n2 := n2;
      dzp.d1 := d1;
      dzp.d2 := d2;
      dzp.k := k;
      dzp.Ts := Ts;
      dzp.method := method;
      dzp.uName := uName;
      dzp.yName := yName;
    end fromFactorization;

  end 'constructor';

encapsulated operator '-'
  function subtract "Subtract two TransferFunctions (zp1 - zp2)"
    import Modelica;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Complex;

    input ZerosAndPoles zp1;
    input ZerosAndPoles zp2;

    protected
    Integer size_z1n1=size(zp1.n1, 1);
    Integer size_z1d1=size(zp1.d1, 1);
    Integer size_z1n2=size(zp1.n2, 1);
    Integer size_z1d2=size(zp1.d2, 1);
    Integer size_z2n1=size(zp2.n1, 1);
    Integer size_z2d1=size(zp2.d1, 1);
    Integer size_z2n2=size(zp2.n2, 1);
    Integer size_z2d2=size(zp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output ZerosAndPoles result "= zp1/zp2";

  algorithm
    if zp1 == zp2 then
      result := ZerosAndPoles(0);
    else
      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,zp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          zp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,zp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          zp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,zp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          zp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,zp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          zp1.d2[i, :]));
      end for;

      p3 := zp1.k*p1 - zp2.k*p2;
      k := 0;
      for i in size(p3.c, 1):-1:1 loop
        if abs(p3.c[i]) > Modelica.Constants.eps then
          k := p3.c[i];
        end if;
      end for;
      numZeros := Polynomial.roots(p3);
      result := ZerosAndPoles(
        numZeros,
        dummy,
        k);

      result.d1 := cat(
        1,
        zp1.d1,
        zp2.d1);
      result.d2 := cat(
        1,
        zp1.d2,
        zp2.d2);

    end if;
  end subtract;

  function negate "Unary minus (multiply transfer function by -1)"
      import ZerosAndPoles =
        Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

     input ZerosAndPoles zp;
     output ZerosAndPoles result(n1=zp.n1, n2=zp.n2, d1=zp.d1, d2=zp.d2, k=-zp.k) "= -zp";
  algorithm
  end negate;
end '-';

  encapsulated operator function '+'
    "Addition of to tarnsfwer functions zp1 + zp2, i.e. parallel connection of two transfer functions (= inputs are the same, outputs of the two systems are added)"

    import Modelica;
    import ZerosAndPoles = Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Complex;

    input ZerosAndPoles zp1;
    input ZerosAndPoles zp2;

  protected
    Integer size_z1n1=size(zp1.n1, 1);
    Integer size_z1d1=size(zp1.d1, 1);
    Integer size_z1n2=size(zp1.n2, 1);
    Integer size_z1d2=size(zp1.d2, 1);
    Integer size_z2n1=size(zp2.n1, 1);
    Integer size_z2d1=size(zp2.d1, 1);
    Integer size_z2n2=size(zp2.n2, 1);
    Integer size_z2d2=size(zp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output ZerosAndPoles result "= zp1+zp2";

  algorithm
    if zp1 == -zp2 then
      result := ZerosAndPoles(0);
    else

      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,zp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          zp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,zp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          zp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,zp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          zp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,zp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          zp1.d2[i, :]));
      end for;

      p3 := zp1.k*p1 + zp2.k*p2;
      k := p3.c[1];
      numZeros := Polynomial.roots(p3);
      result := ZerosAndPoles(
        numZeros,
        dummy,
        k);

      result.d1 := cat(
        1,
        zp1.d1,
        zp2.d1);
      result.d2 := cat(
        1,
        zp1.d2,
        zp2.d2);
    end if;

  end '+';

  encapsulated operator function '*'
    "Multiply two ZerosAndPoles transfer functions (zp1 * zp2)"

    import Modelica;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    input ZerosAndPoles zp1;
    input ZerosAndPoles zp2;

    output ZerosAndPoles result "= zp1 * zp2";
  algorithm
    if zp1 == ZerosAndPoles(0) or zp2 == ZerosAndPoles(0) then
      result := ZerosAndPoles(0);
    else
      result.n1 := cat(
        1,
        zp1.n1,
        zp2.n1);
      result.n2 := cat(
        1,
        zp1.n2,
        zp2.n2);
      result.d1 := cat(
        1,
        zp1.d1,
        zp2.d1);
      result.d2 := cat(
        1,
        zp1.d2,
        zp2.d2);
      result.k := zp1.k*zp2.k;
    end if;

  end '*';

  encapsulated operator function '/'
    "Divide two transfer functions (zp1 / zp2)"

    import Modelica;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    input ZerosAndPoles zp1;
    input ZerosAndPoles zp2;

    output ZerosAndPoles result "= zp1/zp2";

  algorithm
    assert(abs(zp2.k)>100*Modelica.Constants.small,"zp2 in operator \"Modelica_LinearSystems2.TransferFunction.'/'()\" may not be zero");
    if zp1==ZerosAndPoles(0) then
      result := ZerosAndPoles(0);
    else
      result.n1 := cat(1,zp1.n1, zp2.d1);
      result.n2 := cat(1,zp1.n2, zp2.d2);
      result.d1 := cat(1,zp1.d1, zp2.n1);
      result.d2 := cat(1,zp1.d2, zp2.n2);
      result.k := zp1.k/zp2.k;
  end if;

  end '/';

  encapsulated operator function '^' "Integer power of TransferFunction (zp^k)"

    import Modelica;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    input ZerosAndPoles zp;
    input Integer k;

    output ZerosAndPoles result(
      redeclare Real n1[k*size(zp.n1, 1)],
      redeclare Real n2[k*size(zp.n2, 1),2],
      redeclare Real d1[k*size(zp.d1, 1)],
      redeclare Real d2[k*size(zp.d2, 1),2]) "= zp^k";
  protected
    Integer size_n1=size(zp.n1, 1);
    Integer size_d1=size(zp.d1, 1);
    Integer size_n2=size(zp.n2, 1);
    Integer size_d2=size(zp.d2, 1);
  algorithm
    for i in 1:k loop
      result.n1[(i - 1)*size_n1 + 1:i*size_n1] := zp.n1;
      result.d1[(i - 1)*size_d1 + 1:i*size_d1] := zp.d1;
      result.n2[(i - 1)*size_n2 + 1:i*size_n2, :] := zp.n2;
      result.d2[(i - 1)*size_d2 + 1:i*size_d2, :] := zp.d2;
    end for;
    result.k := zp.k^k;

  end '^';

encapsulated operator function '=='
    "Check whether two transfer functions are identical"
    import Modelica;
    import Modelica.Math;
    import Modelica_LinearSystems2.Math.Polynomial;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    input ZerosAndPoles zp1;
    input ZerosAndPoles zp2;
    input Real eps(min=0) = 0
      "Two numbers n1 and n2 are identical if abs(n1-n2) <= eps";

    output Boolean result "= zp1 == zp2";
algorithm
    result := Math.Vectors.isEqual(zp1.n1,zp2.n1,eps) and Math.Vectors.isEqual(zp1.d1,zp2.d1,eps) and Math.Matrices.isEqual(zp1.n2,zp2.n2,eps) and Math.Matrices.isEqual(zp1.d2,zp2.d2,eps) and (zp1.k==zp2.k);

end '==';

  encapsulated operator function 'String'
    "Transform ZerosAndPoles transfer function into a String representation"
    import Modelica;
    import ZerosAndPoles =
      Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    //import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles.Internal;
    import Modelica_LinearSystems2.ZerosAndPoles.Internal;

      input ZerosAndPoles zp
      "ZerosAndPoles transfer function to be transformed in a String representation";
      input Integer significantDigits=6
      "Number of significant digits that are shown";
      input String name="p" "Independent variable name used for printing";
      output String s="";
  protected
      Integer num_order=size(zp.n1, 1) + 2*size(zp.n2, 1);
      Integer den_order=size(zp.d1, 1) + 2*size(zp.d2, 1);
  algorithm
      if num_order == 0 and den_order == 0 then
        s := String(zp.k);

      else
         // construct string for multiplicative factor
        if zp.k <> 1.0 or zp.k == 1.0 and num_order == 0 then
          s := String(zp.k);
          if num_order <> 0 then
            s := s + "*";
          end if;
        end if;

        if num_order <> 0 then
            // construct numerator string
          s := s + Internal.firstOrderToString(
                zp.n1,
                significantDigits,
                name);
          if size(zp.n2, 1) <> 0 then
            s := if size(zp.n1, 1) > 0 then s + "*" +
              Internal.secondOrderToString(
                  zp.n2,
                  significantDigits,
                  name) else s + Internal.secondOrderToString(
                  zp.n2,
                  significantDigits,
                  name);
          end if;
        end if;

        if den_order <> 0 then
            // construct denominator string
          s := s + "/";
          if den_order > 1 then
            s := s + "( ";
          end if;
          s := s + Internal.firstOrderToString(
                zp.d1,
                significantDigits,
                name);
          if size(zp.d2, 1) <> 0 then
            if size(zp.d1, 1) > 0 then
              s := s + "*";
            end if;
            s := s + Internal.secondOrderToString(
                  zp.d2,
                  significantDigits,
                  name);
          end if;
          if den_order > 1 then
            s := s + " )";
          end if;
        end if;
      end if;
    //    end toString;
  end 'String';

  encapsulated function q "Generate the transfer function p"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[1],
      redeclare Real n2[0,2],
      redeclare Real d1[0],
      redeclare Real d2[0,2]);
  algorithm
    dzp.n1[1] := 0;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
q = DiscreteZerosAndPoles.<b>q</b>()
</pre></blockquote>

<h4>Description</h4>
<p>
Generate the complex Laplace variable q=rxp(s*T) as a DiscreteZerosAndPoles transfer function. It can be used for generating like
</p>
<blockquote><pre>
DiscreteZerosAndPoles dzp = q/(q^2 + q + 1)/(q + 1)
</pre></blockquote>
</html>"));
  end q;

  encapsulated package Analysis

  encapsulated function denominatorDegree "Return denominator degree"
      import Modelica;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Integer result;
  algorithm
    result := size(dzp.d1, 1) + 2*size(dzp.d2, 1);
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> DiscreteZerosAndPoles.Analysis.<b>denominatorDegree</b>(dzp)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial constituted by the first and second order polynomials of the DiscreteZeroAndPoles numerator.
</p>

<h4>Example</h4>
<blockquote><pre>
   DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
   Modelica_LinearSystems2.DiscreteZerosAndPoles zp=(p+1)/(p^2+p+1);

   Real dDegree;

<b>algorithm</b>
  dDegree := DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp);
//  dDegree = 2
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.numeratorDegree\">DiscreteZerosAndPoles.Analysis.numeratorDegree</a>
</p>
</html>"));
  end denominatorDegree;

  end Analysis;

  encapsulated package Design
  end Design;

  encapsulated package Plot
  end Plot;

  encapsulated package Conversion

  function toDiscreteTransferFunction
      "Generate a DiscreteTransferFunction object from a DiscreteZerosAndPoles object"

    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Internal;
    import Complex;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output DiscreteTransferFunction dtf;

    protected
    ZerosAndPoles zp=ZerosAndPoles(k=dzp.k, n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2);
    Real k;
    Complex z[:];
    Complex p[:];
    Polynomial pn;
    Polynomial pd;
  algorithm
    (z,p,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
    pn := Polynomial(z)*Polynomial(k);
    pd := Polynomial(p);
    dtf.n := pn.c;
    dtf.d := pd.c;
    dtf.Ts := dzp.Ts;
    dtf.method := dzp.method;
    dtf.uName := dzp.uName;
    dtf.yName := dzp.yName;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> ZerosAndPoles.Conversion.<b>toTransferFunction</b>(zp)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a TransferFunction record
 <blockquote><pre>
           n(s)     b0 + b1*s + ... + bn*s^n
   tf = -------- = --------------------------
           d(s)     a0 + a1*s + ... + an*s^n
 </pre></blockquote>
from a ZerosAndPoles record representated by first and second order numerator and denominator polynomials. The poles and zeros and the gain <tt>k</tt> are computed (<a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs in the TransferFunction constructor.


<h4>Example</h4>
<blockquote><pre>
   ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
   Modelica_LinearSystems2.ZerosAndPoles tf = 1/(p + 3)/(p + 1)


<b>algorithm</b>
  zp:=Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction(tf);
//  zp = 1/( (p + 1)*(p + 2) )
</pre></blockquote>



</html>"));
  end toDiscreteTransferFunction;

  function toStateSpace
      "Transform a ZerosAndPoles object into a StateSpace object"
   //encapsulated function fromZerosAndPoles
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Math.Vectors;
      import Complex = Modelica_LinearSystems2.Math.ComplexAdvanced;
      import Modelica_LinearSystems2.StateSpace;
      import
        Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles.Internal;

    input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
    output StateSpace ss(
      redeclare Real A[ZerosAndPoles.Analysis.denominatorDegree(zp),
        ZerosAndPoles.Analysis.denominatorDegree(zp)],
      redeclare Real B[ZerosAndPoles.Analysis.denominatorDegree(zp),1],
      redeclare Real C[1,ZerosAndPoles.Analysis.denominatorDegree(zp)],
      redeclare Real D[1,1]) "Transfer function in StateSpace SISO form";

    protected
    Real A[2,2] "system matrix of partial 2nd order system";
    Real B[2,1] "input matrix of partial 2nd order system";
    Real C[1,2] "output matrix of partial 2nd order system";
    Real D[1,1] "feedthrough matrix of partial 2nd order system";
    Real a "system 'matrix' of partial 1st order system";
    Real b "input 'matrix' of partial 1st order system";
    Real c "output 'matrix' of partial 1st order system";
    Real d "feedthrough 'matrix' of partial 1st order system";
    Integer nx=max(ZerosAndPoles.Analysis.numeratorDegree(zp),ZerosAndPoles.Analysis.denominatorDegree(zp));
    Integer n_num1=size(zp.n1, 1);
    Integer n_num2=size(zp.n2, 1);
    Integer n_den1=size(zp.d1, 1);
    Integer n_den2=size(zp.d2, 1);
    Integer n_num=n_num1 + 2*n_num2;
    Integer n_den=n_den1 + 2*n_den2;

    Integer i_d=if n_num2 > n_den2 then 2*(n_num2 - n_den2) + 1 else 1;
    Integer i_k=if n_num2 > n_den2 then n_den2 - (n_num2 - n_den2) else n_den2;
    Integer i;
    Integer ili;
    Real num[nx,2]=[zp.n2; [zp.n1,zeros(n_num1)]; zeros(max(0,nx - n_num2 - n_num1), 2)]
        "Numerator matrix, in order that indices are defined in all situations in all if clauses";
    Real den[nx,2]=[zp.d2; [zp.d1,zeros(n_den1)]; zeros(max(0,nx - n_den2 - n_den1), 2)]
        "Denominator matrix, in order that indices are defined in all situations in all if clauses";
    Real k[i_k + n_den1](each fixed=false)
        "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
    Real k_total;

    Boolean dZero=true;

  algorithm
    assert(n_num <= n_den,
      "ZerosAndPoles transfer function is not proper as required from StateSpace system:\n"
       + "  numerator degree (= " + String(n_num) +
      ") <= denominator degree (= " + String(n_den) + ") required.");

    if n_den > 0 then

      for i in 1:max(n_den2, n_num2) loop
        // State space systems of order 2
        if i <= n_den2 then
          if i <= n_num2 then
              // State space system in form (1)
            k[i] := Internal.scaleFactor2(
                num[i, 1],
                num[i, 2],
                den[i, 1],
                den[i, 2]);
          elseif i - n_num2 + 1 <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
            k[i] := Internal.scaleFactor2(
                num[i, 1] + num[i + 1, 1],
                num[i, 1]*num[i + 1, 1],
                den[i, 1],
                den[i, 2]);
          elseif i - n_num2 == n_num1 then
              // State space system in form (2) with 1 first order numerator polynomial
            k[i] := Internal.scaleFactor2(
                0,
                num[i, 1],
                den[i, 1],
                den[i, 2]);
          else
              // State space system in form (3)
            k[i] := Internal.scaleFactor2(
                0,
                0,
                den[i, 1],
                den[i, 2]);
          end if;
        else
           // State space system in form (1) with 2 first order denominator polynomials
          k[i] := Internal.scaleFactor2(
              num[i, 1],
              num[i, 2],
              den[i, 1] + den[i + 1, 1],
              den[i, 1]*den[i + 1, 1]);
        end if;
      end for;

      for i in i_d:n_den1 loop
        // State space systems of order 1
        if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
           // State space system in form (4)
          k[i_k + i] := Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
            n_num2) + i), 1], den[n_den2 + i, 1]);
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
           // State space system in form (4)
          k[i_k + i] := Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
            1], den[n_den2 + i, 1]);
        else
           // State space system in form (5)
          k[i_k + i] := Internal.scaleFactor1(0, den[n_den2 + i, 1]);
        end if;
      end for;

      k_total := zp.k/product(k);

      ss.A := zeros(nx, nx);
      ss.B := zeros(nx, 1);
      ss.C := zeros(1, nx);
      ss.D := zeros(1, 1);

   // Calculation of matrices A, B, C, D
   //first elements of A, B, C and D

      if max(n_den2, n_num2) > 0 then
        ili := i_d;
        A[1, :] := {0,1};
        B[1, 1] := 0;
            // Construct state space systems of order 2
        if 1 <= n_den2 then
          A[2, :] := {-den[1, 2],-den[1, 1]};
          B[2, 1] := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then den[1, 2] else 1;
          if 1 <= n_num2 then
                 // State space system in form (1)
            C := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;

          elseif 1 - n_num2 + 1 <= n_num1 then
                 // State space system in form (1) with 2 first order numerator polynomials

            C := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then  k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;
          elseif 1 - n_num2 == n_num1 then
                 // State space system in form (2) with 1 first order numerator polynomial

            C := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then k[1]*[num[1, 1],1]/den[1, 2] else k[1]*[num[1, 1],1];
            D := [0];
            dZero := dZero and true;
          else
                 // State space system in form (3)

            C := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then  k[1]*[1,0]/den[1, 2] else k[1]*[1,0];
            D := [0];
            dZero := dZero and true;
          end if;
        else
              // State space system in form (1) with 2 first order denominator polynomials

          A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
          B[2, 1] := den[1, 1]*den[1, 1];

          C := k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/den[1, 1]
            /den[1, 1];
          D := [k[1]];
          dZero := false;
        end if;
        ss.A[1:2, 1:2] := A;
        ss.B[1:2, 1] := vector(B);
        ss.C[1, 1:2] := vector(C);
        ss.D := D;

      else
        ili := max(2, i_d);
     // Construct state space systems of order 1
        a := -den[1, 1];
        b := if abs(den[1, 1]+1)>Modelica.Constants.eps then den[1,1]+1 else if n_num1>0 then num[1,1] else 1;

        if 1 <= n_num1 then
              // State space system in form (4)
          c := if abs(den[1, 1]+1)>Modelica.Constants.eps then k[1]*(num[1, 1] - den[1, 1])/(den[1, 1]+1) else k[1];
          d := k[1];
          dZero := false;
        else
              // State space system in form (5)
          c := if abs(den[1, 1]+1)>Modelica.Constants.small then k[1]/(den[1, 1]+1) else k[1];
          d := 0;
          dZero := dZero and true;
        end if;
        ss.A[1, 1] := a;
        ss.B[1, 1] := b;
        ss.C[1, 1] := c;
        ss.D[1, 1] := d;

      end if;

   /// for i=2 to degree(system)
      A[1, :] := {0,1};
      B[1, 1] := 0;

      for i in 2:max(n_den2, n_num2) loop
           // Construct state space systems of order 2
        if i <= n_den2 then
          A[2, :] := {-den[i, 2],-den[i, 1]};
          B[2, 1] := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then den[i, 2] else 1;

          if i <= n_num2 then
                 // State space system in form (1)

            C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif i - n_num2 + 1 <= n_num1 then
                 // State space system in form (1) with 2 first order numerator polynomials
            C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then k[i]*[num[i, 1]*num[i + 1, 1] - den[i, 2],num[i, 1] + num[i + 1, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[i, 1]*num[i + 1, 1] - den[i, 2],num[i, 1] + num[i + 1, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif i - n_num2 == n_num1 then
                 // State space system in form (2) with 1 first order numerator polynomial
            C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then k[i]*[num[i, 1],1]/den[i, 2] else k[i]*[num[i, 1],1];
            D := [0];
            dZero := dZero and true;

          else
                 // State space system in form (3)
            C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then k[i]*[1,0]/den[i, 2] else k[i]*[1,0]/den[i, 2];
            D := [0];
            dZero := dZero and true;

          end if;

        else
              // State space system in form (1) with 2 first order denominator polynomials
          A[2, :] := {-(den[i, 1]*den[i + 1, 1]),-(den[i, 1] + den[i + 1, 1])};
          B[2, 1] := den[i, 1]*den[i, 1];

          C := k[i]*[num[i, 2] - (den[i, 1]*den[i + 1, 1]),num[i, 1] - (den[i, 1] + den[i + 1, 1])]/
            den[i, 1]/den[i, 1];
          D := [k[i]];
          dZero := false;
        end if;

        ss.A[2*i, 1:2*i - 2] := B[2, 1]*ss.C[1, 1:2*i - 2];
        ss.A[2*i - 1:2*i, 2*i - 1:2*i] := A;
        ss.B[2*i, 1] := if dZero then 0 else B[2, 1]*ss.D[1, 1];
        ss.C[1, 1:2*i - 2] := if dZero then fill(0, 2*i - 2) else D[1, 1]*ss.C[
          1, 1:2*i - 2];
        ss.C[1, 2*i - 1:2*i] := vector(C);
        ss.D := D*ss.D;
      end for;

   //  for i in max(2,i_d):n_den1 loop

      for i in ili:n_den1 loop
           // Construct state space systems of order 1
        a := -den[n_den2 + i, 1];
        b := if abs(den[n_den2 + i, 1]+1)>Modelica.Constants.eps then den[n_den2 + i, 1]+1 else 1.0;
        if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
              // State space system in form (4)

          c := if abs(den[n_den2 + i, 1]+1)>Modelica.Constants.eps then k[i_k + i]*(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/(den[n_den2 + i, 1]+1) else 1.0;

          d := k[i_k + i];

          dZero := false;
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
              // State space system in form (4)

          c := if abs(den[n_den2 + i, 1]+1)>Modelica.Constants.eps then k[i_k + i]*(num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/(den[n_den2 + i, 1]+1) else 1.0;
          d := k[i_k + i];
          dZero := false;
        else

              // State space system in form (5)

          c := if abs(den[n_den2 + i, 1]+1)>Modelica.Constants.eps then k[i_k + i]/(den[n_den2 + i, 1]+1) else k[i_k + i];
          d := 0;
          dZero := dZero and true;
        end if;

        ss.A[2*n_den2 + i, 1:2*n_den2 + i - 1] := b*ss.C[1, 1:2*n_den2 + i - 1];
        ss.A[2*n_den2 + i, 2*n_den2 + i] := a;
        ss.B[2*n_den2 + i, 1] := if dZero then 0 else b*ss.D[1, 1];
        ss.C[1, 1:2*n_den2 + i - 1] := if dZero then fill(0, 2*n_den2 + i - 1) else
                d*ss.C[1, 1:2*n_den2 + i - 1];
        ss.C[1, 2*n_den2 + i] := c;
        ss.D := if dZero then [0] else d*ss.D;

      end for;
      ss.C := k_total*ss.C;
      ss.D := k_total*ss.D;
    else
      ss := Modelica_LinearSystems2.StateSpace(zp.k);
    end if;

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = ZerosAndPoles.Conversion.<b>toStateSpace</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
This function transforms a zeros-poles-gain system representation into state space representation.
To achieve well numerical condition the ZerosAndPoles transfer function is transformed into state space
form by creating first and second order blocks that are connected
together in series. Every block is represented in controller
canonical form and scaled such that the gain from the input
of this block to its output is one (i.e. y(p=0) = u(p=0)),
if this is possible. Details are given below.
</p>

<h4>Algorithmic details</h4>
<p>
The ZerosAndPoles transfer function is defined as:
</p>
<blockquote><pre>
         product(p + n1[i]) * product(p^2 + n2[i,1]*p + n2[i,2])
y = k * --------------------------------------------------------- * u
         product(p + d1[i]) * product(p^2 + d2[i,1]*p + d2[i,2])
</pre></blockquote>
<p>
This is treated as a series connection of first and second order
systems. If size(n1) == size(d1) and size(n2) == size(d2)
this gives the following sequence of operations:
</p>
<blockquote><pre>
        p^2 + n2[1,1]*p + n2[1,2]
  y_1 = ------------------------- * u
        p^2 + d2[1,1]*p + d2[1,2]
&nbsp;
        p^2 + n2[2,1]*p + n2[2,2]
  y_2 = ------------------------- * y_1
        p^2 + d2[2,1]*p + d2[2,2]
&nbsp;
     ...
&nbsp;
        p + n1[..]
  y_n = ---------- * y_(n-1)
        p + d1[..]
&nbsp;
    y = k*y_n
</pre></blockquote>
<p>
Based on this representation, evrey block with transfer function G(p) could be transformed into
</p>
<blockquote><pre>
G(p) = k * F(p)
</pre></blockquote>
<p>
with F(p) has unit gain. This leads to representations of the forms
</p>
<blockquote><pre>
         a2 + a1*p + p^2       a2      b2 + a1*b2/a2*p + b2/a2*p^2
G(p) = -------------------- = ---- * ------------------------------ = k * F(p),  k = a2/b2  (1)
         b2 + b1*p + p^2       b2           b2 + b1*p + p^2
</pre></blockquote>
<p>
for second order systems and
</p>
<blockquote><pre>
         a + p     a     b + b/a*p
G(p) = -------- = --- * ---------- = k * F(p),   k = a/b
         b + p     b      b + p
</pre></blockquote>
<p>
for first order systems respectively.
</p>

<p>
The complete system is now considered as the series connections of all the single unit gain transfer functions and an overall gain k with
</p>
<blockquote><pre>
k = product(ki).
</pre></blockquote>
<p>
In the general case, the following system structures
and the corresponding state space systems can appear
(note, 'c' is the reciprocal local gain 1/k):
</p>
<blockquote><pre>
(1)
          a2 + a1*p + p^2           der(x1) = x2
    y = ---------------------  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
          b2 + b1*p + p^2                 y = c*((a2-b2)/b2*x1 + (a1-b1)/b2*x2 + u),  c = b2/a2
&nbsp;
(2)
             p + a                 der(x1) = x2
    y = ---------------- * u  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                  y = k*(a/b2*x1 +x2/b2),  c = b2/a
&nbsp;
(3)
               1                  der(x1) = x2
    y = --------------- *u   -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                 y = c*x1/b2,  c = b2
&nbsp;
(4)
       a + p                       der(x) = -b*x + b*u
   y = ----- * u             -->        y = c*((a-b)/b*x + u),  c = b/a
       b + p
&nbsp;
(5)
         1
   y = ----- * u             -->   der(x) = -b*x + b*u
       b + p                            y = x,  c = b

</pre></blockquote>
<p>
If the sizes of the numerator and denominator polynomials
do not match, the small systems are built in the
following way:
</p>
<blockquote><pre>
(1) Build systems of form (1) by combining
    - 1 d2 and 1 n2
      (= 1 second order denominator and 1 second order numerator) or
    - 1 d2 and 2 n1 or
    - 2 d1 and 1 n2
(2) Build at most one system of form (2) by combining
    - 1 d2 and 1 n2
(3) Build systems of form (3) by
    - 1 d2
(4) Build systems of form (4) by combining
    - 1 d1 and 1 n1
(5) Build systems of form (5) by
    - 1 d1
</pre></blockquote>
<p>
The numeric properties of the resulting state space system
depends on which first and second order polynomials are
combined and connected together. From a numerical point of view, it
would therefore be useful to combine the polynomials
based on the numeric values of the polynomial coefficients,
(e.g., in a first step the polynomials could be sorted
according to their cut-off frequency).
</p>
<p>
However, this has the disadvantage that the structure of the
resulting state space system depends on the numeric
values of the polynomial coefficients. Since Modelica
environments perform symbolic pre-processing on equations,
this would mean that a change of a polynomial coefficient
requires to newly compile the state space system.
</p>
<p>
If, on the other hand, the structure of the state
space system depends only on dimension information
of the n1,n2,d1,d2 arrays, then the polynomial coefficients
can be changed without a new translation of the model.
This is the major reason why the structure of the
state space system in the implementation of this block
is based only on dimension information.
</p>
<p>
This is, e.g., not critical for the provided filters:
The dimension of the n1,n2,d1,d2 arrays depend for
filters only on the filter characteristics
(Bessel, Butterworth etc.), the filter type (low pass,
high pass etc.) and on the filter order. If any
of this data is changed, the model has to be
newly compiled. All the other filter data, such as
cut-off frequency or ripple amplitude, can be changed
without re-compilation of the model.
The ZerosAndPoles transfer function is now constructed
for the filters in such a way that the filter zeros
and poles are appropriately sorted to give better
numerical properties.
</p>
<p>
Another alternative implementation of the state
space system would be to use the function controller canonical
form that directly results from the transfer function.
The severe disadvantage
of this approach is that the structure of the state
space system from above is lost for the symbolic preprocessing.
If, e.g., index reduction has to be applied (e.g. since a
filter is used to realize a non-linear inverse model),
then the tool cannot perform the index reduction.
Example:
</p>
<p>
Assume, a generic first order state space system
is present
</p>
<blockquote><pre>
<b>der</b>(x) = a*x + b*u
     y = c*x + d*u
</pre></blockquote>
<p>
and the values of the scalars a,b,c,d are parameters
that might be changed before the simulation starts.
If y has to be differentiated symbolically during code
generation, then
</p>
<blockquote><pre>
<b>der</b>(y) = c*<b>der</b>(x) + d*<b>der</b>(u)
<b>der</b>(x) = a*x + b*u
</pre></blockquote>
<p>
As a result, u needs to be differentiated too, and this
might not be possible and therefore translation might fail.
</p>
<p>
On the other hand, if the first order system is
defined to be a low pass filter and the state space
system is generated by keeping this structure, we have
(see form (5) above):
</p>
<blockquote><pre>
<b>der</b>(x) = -b*x + u
      y = x
</pre></blockquote>
<p>
Differentiating y symbolically leads to:
</p>
<blockquote><pre>
<b>der</b>(y) = <b>der</b>(x)
<b>der</b>(x) = -b*x + u
</pre></blockquote>
<p>
Therefore, in this case, the derivative of u is not
needed and the tool can continue with the symbolic
processing.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2 + p +1);

<b>algorithm</b>
  ss := Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace(zp);
// ss.A = [0, 1; -1, -1],
// ss.B = [0; 1],
// ss.C = [1, 1],
// ss.D = [0],
</pre></blockquote>

</html>"));
  end toStateSpace;

  end Conversion;

  encapsulated package Import

  function fromModel
      "Generate a ZerosAndPoles record array from a state space representation resulted from linearization of a model"

    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2.StateSpace;
    import ZerosAndPoles =
        Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;
    import Simulator = DymolaCommands.SimulatorAPI;

    input String modelName "Name of the Modelica model";
    input Real T_linearize=0
        "Point in time of simulation to linearize the model";
    input String fileName="dslin" "Name of the result file";

    protected
    String fileName2=fileName + ".mat";
    Boolean OK1=Simulator.simulateModel(
          problem=modelName,
          startTime=0,
          stopTime=T_linearize);
    Boolean OK2=Simulator.importInitial("dsfinal.txt");
    Boolean OK3=Simulator.linearizeModel(
          problem=modelName,
          resultFile=fileName,
          startTime=T_linearize,
          stopTime=T_linearize + 1);
    Integer xuy[3] = ReadSystemDimension(fileName2, "ABCD");
    Integer nx = xuy[1];
    Integer nu = xuy[2];
    Integer ny = xuy[3];
    Real ABCD[nx + ny,nx + nu]=Streams.readRealMatrix(
          fileName2,
          "ABCD",
          nx + ny,
          nx + nu);
    String xuyName[nx + nu + ny]=readStringMatrix(
          fileName2,
          "xuyName",
          nx + nu + ny);

    StateSpace result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "= model linearized at initial point";
    public
    output ZerosAndPoles zp[:,:];
  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    result.uNames := xuyName[nx + 1:nx + nu];
    result.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    result.xNames := xuyName[1:nx];

    zp := StateSpace.Conversion.toZerosAndPolesMIMO(result);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> ZerosAndPoles.Import.<b>fromModel</b>(modelName, T_linearize, fileName)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Generate a matrix of ZerosAndPoles data records by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation. The system is genrated by using <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Import.fromFile\">StateSpace.Import.fromFile</a> followed by a conversion from sate space to transfer function representation.

<h4>Example</h4>
<blockquote><pre>
   String modelName = \"Modelica_LinearSystems2.Examples.DoublePendulum\";
   Real T_linearize = 5;


<b>algorithm</b>
  zp = Modelica_LinearSystems2.ZerosAndPoles.Import.fromModel(modelName, T_linearize);

//  zp =[0.157605*(p + 0.706559)*(p + 12.3798)*(p^2-7.34273*p + 18.674)/( (p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572)*(p^2 + 2.07022e-015*p + 3.38074e-015) );
         0.157605*(p + 0.706559)*(p + 12.3798)*(p^2-7.34273*p + 18.674)/( (p-1.94349e-015)*(p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572) );
        -0.166305*(p^2-1.20297*p + 3.48327)/( (p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572) );
        -0.166305*p*(p^2-1.20297*p + 3.48327)/( (p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572) );
         0.283325*(p-5.23615)*(p + 0.551929)/( (p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572) );
         0.283325*p*(p-5.23615)*(p + 0.551929)/( (p + 0.829834)*(p + 10.6304)*(p^2-7.27298*p + 18.1572) )]

</pre></blockquote>
</html>"));
  end fromModel;

    encapsulated function fromFile
      "Generate a ZerosAndPoles record by reading the polynomial coefficients or zeros and poles from a file"
      //import ZerosAndPoles = Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2;
      import Complex = Modelica_LinearSystems2.Math.ComplexAdvanced;
      import Modelica;
      import Modelica_LinearSystems2.DataDir;

      input String fileName=DataDir + "zp.mat"
        "Name of the zeros and poles data file"
        annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));

    protected
      Integer n1n2d1d2[4]=if ZerosAndPoles.Internal.checkRepresentation(
          fileName) then ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(
          fileName) else ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(
          fileName);
      Integer n1=n1n2d1d2[1];
      Integer n2=n1n2d1d2[2];
      Integer d1=n1n2d1d2[3];
      Integer d2=n1n2d1d2[4];
      Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
      Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
      output ZerosAndPoles zp(
        n1=fill(0, n1),
        n2=fill(
              0,
              n2,
              2),
        d1=fill(0, d1),
        d2=fill(
              0,
              d2,
              2));
    algorithm
    //Whenever this function becomes operational the code must be rewritten if fromFile_pc2 and fromFile_zp2 are in the 'constructor'

      zp := if ZerosAndPoles.Internal.checkRepresentation(fileName) then ZerosAndPoles.Internal.fromFile_zp( fileName) else ZerosAndPoles.Internal.fromFile_pc(
        fileName);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> ZerosAndPoles.Import.<b>fromFile</b>(fileName)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Reads and loads a zeros-and-poles transfer function from a mat-file <tt>fileName</tt>. The file must contain either the set of variables n1, n2, d1, d2, and k with the associated first and second order polynomials or the variables p, z, and k with the poles and zeros, written in two column arrays with real and imaginary in the first and second column respectively. The variable k is the real gail in both cases.


<h4>Example</h4>
<blockquote><pre>


<b>algorithm</b>
  zp:=Modelica_LinearSystems2.ZerosAndPoles.Import.fromFile(\"zp.mat\", \"n\", \"d\");
//  zp = (p^2 + 2*p + 3)/(p + 2)/(p^2 + 2*p + 2)
</pre></blockquote>


</html>"));
    end fromFile;

  end Import;

  encapsulated package Internal
    "Internal library of record Filter (should not be directly used by user)"
    extends Modelica.Icons.InternalPackage;

    import Modelica;
    import Modelica_LinearSystems2;

    function numberOfRealZeros2 "Calculate number of real zeros"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input DiscreteTransferFunction dtf "DiscreteTransferFunction";
    output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(dtf.n)));
    algorithm
    end numberOfRealZeros2;

    function numberOfRealPoles "Calculate number of real poles"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input DiscreteTransferFunction dtf "TransferFunction";
      output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(
          dtf.d)));
    algorithm
    end numberOfRealPoles;

  encapsulated function scaleFactor1
      "Return scale factor for first order block"
      import Modelica;
    input Real n "(z+n)/(z+d)";
    input Real d "(z+n)/(z+d)";
    input Real small=100*Modelica.Constants.eps;
    output Real k "= (1+d)/(1+n), if d,n are not zero, otherwise special cases";
  algorithm
    k := if abs(d+1) > small and abs(d+1) > small then abs(d+1)/abs(n+1) else 1;
  end scaleFactor1;

    function scaleFactor2 "Return scale factor for second order block"
      import Modelica;
    input Real n1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real n2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real small=100*Modelica.Constants.eps;
    output Real k "= d2/n2, if d2,n2 are not zero, otherwise special cases";
    algorithm
    k := (if abs(d2) > small then abs(d2) else (if abs(d1) > small then abs(
      d1) else 1))/(if abs(n2) > small then abs(n2) else (if abs(n1) > small then
            abs(n1) else 1));
    end scaleFactor2;

  end Internal;

  annotation (
    defaultComponentName="filter",
    Documentation(info="<html>
<p>
This record defines a transfer function by its zeros, poles and a gain:
</p>
<pre>         product(p - z[i])
  y = k*------------------- * u
         product(p - n[i])
</pre>
<p>
where z[:] is a Complex vector of zeros, n[:] is a Complex
vector of poles and k is an additional multiplicative factor.
The elements of the two Complex vectors must either be real
numbers or conjugate complex pairs (in order that their product
results in a polynomial with Real coefficients).
</p>
<p>
In the record, the zeros and poles are transformed
into a product of first and second order polynomials.
The data structure is especially useful in applications where first and
second order polynomials are naturally occurring, e.g., as
for <b>filters</b>. In fact, via function
<a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Design.filter\">ZerosAndPoles.Design.filter</a>, a
ZeroAndPole transfer function is generated from
<b>low</b> and <b>high pass</b> analog filters
(<b>CriticalDamping</b>, <b>Bessel</b>, <b>Butterworth</b>, <b>Chebyshev</b>).
The filters are available in <b>normalized</b> (default) and non-normalized form.
In the normalized form, the amplitude of the filter transfer function
at the cutoff frequency is 3 dB.
</p>
<p>
A ZeroAndPole transfer function is internally stored by the coefficients
of first and second order polynomials, and by an additional
multiplicative factor k:
</p>
<pre>         product(s + n1[i]) * product(p^2 + n2[i,1]*p + n2[i,2])
  y = k*---------------------------------------------------------
         product(p + d1[i]) * product(p^2 + d2[i,1]*p + d2[i,2])
</pre>
<p>
Note, the degrees of the numerator and denominator
polynomials are given as:
</p>
<pre>
   degree of numerator   = size(n1,1) + 2*size(n2,1);
   degree of denominator = size(d1,1) + 2*size(d2,1);
</pre>
<p>
Example:
</p>
<pre>                          (p+1)
  zp = 4* -------------------------------------
           (p - 1)*(p - (2+j*3))*(p - (2-j*3))
</pre>
<p>
with j=sqrt(-1), is defined as
</p>
<pre>
   <b>import</b> Modelica_LinearSystems2.Math.Complex;
   <b>import</b> Modelica_LinearSystems2.ZerosAndPoles;

   zp = ZerosAndPoles(z = {Complex(-1,0)},
                      p = {Complex(1,0),
                           Complex(2,3),
                           Complex(2,-3)},
                           k=4);
</pre>
</html>"));
end DiscreteZerosAndPoles;
