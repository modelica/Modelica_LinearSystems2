within Modelica_LinearSystems2;
operator record ZerosAndPoles
  "Continuous zeros and poles description of a single input, single output system (data + operations)"

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
  String uName="u" "Name of input signal" annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal" annotation(Dialog(group="Signal names"));

  encapsulated operator 'constructor'
    "Collection of operators to construct a ZerosAndPoles data record"

    import Modelica;
    import Modelica_LinearSystems2;

    encapsulated function fromReal
      "Generate a ZerosAndPoles data record from a real value"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input Real r "Value of real variable";
      input String uName="" "Input name";
      input String yName="" "Output name";
      output ZerosAndPoles zp(
        redeclare Real n1[0],
        redeclare Real n2[0,2],
        redeclare Real d1[0],
        redeclare Real d2[0,2]) "= r";

    algorithm
      zp.k := r;
      zp.uName := uName;
      zp.yName := yName;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = ZerosAndPoles&apos;constructor&apos;.<b>fromReal</b>(r)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles record zp from a real value, i.e. a without dynamics:
</p>
<blockquote><pre>
y = r*u
</pre></blockquote>
<p>
Therefore, the record is defined by
</p>
<blockquote><pre>
zp.k = r;
zp.n1 = fill(0,1);
zp.n2 = fill(0,1,2);
zp.d1 = fill(0,1);
zp.d2 = fill(0,1,2);
</pre></blockquote>
</html>"));
    end fromReal;

    encapsulated function fromZerosAndPoles
      "Generate a ZerosAndPoles data record from a set of zeros and poles"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Internal;
      import Complex;

      input Complex z[:]=fill(Complex(0,0), 0)
        "Zeros (Complex vector of numerator zeros)";
      input Complex p[:]=fill(Complex(0,0), 0)
        "Poles (Complex vector of denominator zeros)";
      input Real k=1.0 "Constant multiplied with transfer function";
      input String uName="" "Input name";
      input String yName="" "Output name";
      output ZerosAndPoles zp(
        redeclare Real n1[Internal.numberOfRealZeros(z)],
        redeclare Real n2[integer((size(z, 1) - Internal.numberOfRealZeros(z))/2),
        2],
        redeclare Real d1[Internal.numberOfRealZeros(p)],
        redeclare Real d2[integer((size(p, 1) - Internal.numberOfRealZeros(p))/2),
        2]) "ZerosAndPoles transfer functions of the zeros, poles and k";

    protected
      Integer n_n1=size(zp.n1, 1);
      Integer n_d1=size(zp.d1, 1);
      Integer n_n2=size(zp.n2, 1);
      Integer n_d2=size(zp.d2, 1);
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
        zp.n1[i] := -z_reordered[i].re;
      end for;

      j := 1;
      for i in n_n1 + 1:2:size(z, 1) loop
        zp.n2[j, :] := {-2*z_reordered[i].re,z_reordered[i].re^2 + z_reordered[i].im
          ^2};
        j := j + 1;
      end for;

      // Denominator
      for i in 1:n_d1 loop
        zp.d1[i] := -p_reordered[i].re;
      end for;

      j := 1;
      for i in n_d1 + 1:2:size(p, 1) loop
        zp.d2[j, :] := {-2*p_reordered[i].re,p_reordered[i].re^2 + p_reordered[i].im
          ^2};
        j := j + 1;
      end for;

      zp.k := k;
      zp.uName := uName;
      zp.yName := yName;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = ZerosAndPoles&apos;constructor&apos;.<b>fromPolesAndZeros</b>(z, p, k)
   or
zp = ZerosAndPoles&apos;constructor&apos;.<b>fromPolesAndZeros</b>(z, p, k, uName, yName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles transfer function from denominator
and numerator zeros, as well as a gain.
</p>
<p>
Since only transfer functions with real coefficients are supported,
complex roots must be defined as conjugate complex pairs.
It is required that complex conjugate pairs must directly
follow each other as above. An error occurs if this is not the case.
</p>

<h4>Example</h4>
<blockquote><pre>
                       (s+1)
zp = 4 * -------------------------------------
          (s - 1)*(s - (2+j*3))*(s - (2-j*3))
</pre></blockquote>
<p>
with j=sqrt(-1), is defined as
</p>
<blockquote><pre>
  <b>import</b> Modelica_LinearSystems2.Math.Complex;
  <b>import</b> Modelica_LinearSystems2.ZerosAndPoles;

  zp = ZerosAndPoles(z = {Complex(-1,0)},
                     p = {Complex(1,0),
                          Complex(2,3),
                          Complex(2,-3)},
                          k=4);
</pre></blockquote>
</html>"));
    end fromZerosAndPoles;

    function fromTransferFunction =
      Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles
      "Generate a ZerosAndPoles data record from a transfer function"
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
zp = ZerosAndPoles&apos;constructor&apos;.<b>fromTransferFunction</b>(tf)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles record zp from a transfer function tf.
For the simplicity of implementation, this function directly extends from
<a href=\"Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles\">TransferFunction.Conversion.toZerosAndPoles</a>.
</p>
</html>"));

    encapsulated function fromFactorization
      "Generate a ZerosAndPoles data record from first and second order polynomials"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input Real n1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order numerator polynomials"
        annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real n2[:,2]=fill(0, 0, 2)
        "[p,p^0] coefficients of 2nd order numerator polynomials"
        annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order denominator polynomials"
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d2[:,2]=fill(0, 0, 2)
        "[p,p^0] coefficients of 2nd order denominator polynomials"
        annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real k=1.0 "Multiplicative factor of transfer function"
        annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input String uName="" "Input name";
      input String yName="" "Output name";
      output ZerosAndPoles zp(
        redeclare Real n1[size(n1, 1)],
        redeclare Real n2[size(n2, 1),2],
        redeclare Real d1[size(d1, 1)],
        redeclare Real d2[size(d2, 1),2]) "ZerosAndPoles transfer function";
    algorithm
      zp.n1 := n1;
      zp.n2 := n2;
      zp.d1 := d1;
      zp.d2 := d2;
      zp.k := k;
      zp.uName := uName;
      zp.yName := yName;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
zp = ZerosAndPoles&apos;constructor&apos;.<b>fromFactorization</b>(n1, n2, d1, d2, k, uName, yName)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles record zp from first and second order polynomials.
</p>
</html>"));
    end fromFactorization;

    annotation (Documentation(info="<html>
<p>This package contains the default constructors for a data record of zeros-and-poles transfer function.</p>
</html>"), Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25.0),
          Rectangle(
            lineColor={128,128,128},
            fillPattern=FillPattern.None,
            extent={{-100,-100},{100,100}},
            radius=25.0)}));
  end 'constructor';

encapsulated operator '-'
  "Collection of operators for subtraction of zeros and poles descriptions"
  import Modelica;

  function subtract "Subtract two zeros and poles descriptions (zp1 - zp2)"
    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;

    input ZerosAndPoles zp1 "Zeros-and-poles data record 1";
    input ZerosAndPoles zp2 "Zeros-and-poles data record 2";

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
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 + size_z2d2));
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

  function negate "Unary minus (multiply zeros and poles description by -1)"
      import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp "Zeros-and-poles data record";
    output ZerosAndPoles result(n1=zp.n1, n2=zp.n2, d1=zp.d1, d2=zp.d2, k=-zp.k) "= -zp";
  algorithm
  end negate;
    annotation (Documentation(info="<html>
<p>This package contains operators for subtraction of zeros and poles descriptions. </p>
</html>"), Icon(graphics={
        Rectangle(
          lineColor={200,200,200},
          fillColor={248,248,248},
          fillPattern=FillPattern.HorizontalCylinder,
          extent={{-100,-100},{100,100}},
          radius=25.0),
        Line(
          points={{-50,0},{50,0}},
          color={0,0,0},
          smooth=Smooth.None),
        Rectangle(
          lineColor={128,128,128},
          fillPattern=FillPattern.None,
          extent={{-100,-100},{100,100}},
          radius=25.0)}));
end '-';

  encapsulated operator function '+'
    "Addition of two zeros and poles descriptions zp1 + zp2, i.e. parallel connection of two transfer functions (= inputs are the same, outputs of the two systems are added)"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;

    input ZerosAndPoles zp1 "Zeros-and-poles data record 1";
    input ZerosAndPoles zp2 "Zeros-and-poles data record 2";

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
    "Multiply two zeros and poles descriptions (zp1 * zp2)"

    import Modelica;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp1 "Zeros-and-poles data record 1";
    input ZerosAndPoles zp2 "Zeros-and-poles data record 2";

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
    "Divide two zeros and poles descriptions (zp1 / zp2)"
    import Modelica;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp1 "Zeros-and-poles data record 1";
    input ZerosAndPoles zp2 "Zeros-and-poles data record 2";
    output ZerosAndPoles result "Result = zp1/zp2";

  algorithm
    assert(abs(zp2.k)>100*Modelica.Constants.small,"Record zp2 in operator \"Modelica_LinearSystems2.TransferFunction.'/'()\" may not be zero");
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

  encapsulated operator function '^'
    "Integer power of zeros and poles description (zp^k)"

    import Modelica;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp "Zeros-and-poles data record";
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
    "Check whether two zeros and poles descriptions are identical"
  import Modelica;
  import Modelica.Math;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input ZerosAndPoles zp1 "Zeros-and-poles data record 1";
  input ZerosAndPoles zp2 "Zeros-and-poles data record 2";
  input Real eps(min=0) = 0
      "Two numbers n1 and n2 are identical if abs(n1-n2) <= eps";
  output Boolean result "= zp1 == zp2";

algorithm
  result := Math.Vectors.isEqual(zp1.n1,zp2.n1,eps) and Math.Vectors.isEqual(zp1.d1,zp2.d1,eps) and Math.Matrices.isEqual(zp1.n2,zp2.n2,eps) and Math.Matrices.isEqual(zp1.d2,zp2.d2,eps) and (zp1.k==zp2.k);

end '==';

  encapsulated operator function 'String'
    "Transform zeros and poles description into a String representation"
    import Modelica;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp
      "Zeros-and-poles data record to be transformed in a String representation";
    input Integer significantDigits=6
      "Number of significant digits that are shown";
    input String name="p" "Independent variable name used for printing";
      //input Boolean normalized = true;
    output String s="";
  protected
    Boolean normalized = false;
    Real gain=1.0;
    Integer num_order=size(zp.n1, 1) + 2*size(zp.n2, 1);
    Integer den_order=size(zp.d1, 1) + 2*size(zp.d2, 1);
    String sn1;
    String sn2;
    String sd1;
    String sd2;
    Real kn1;
    Real kn2;
    Real kd1;
    Real kd2;
  algorithm
    if num_order == 0 and den_order == 0 then
      s := String(zp.k);
      return;
    end if;

    // construct numerator and denominator strings
    (sn1,kn1) :=ZerosAndPoles.Internal.firstOrderToString( zp.n1, significantDigits, name, normalized);
    (sn2,kn2) :=ZerosAndPoles.Internal.secondOrderToString(zp.n2, significantDigits, name, normalized);
    (sd1,kd1) :=ZerosAndPoles.Internal.firstOrderToString( zp.d1, significantDigits, name, normalized);
    (sd2,kd2) :=ZerosAndPoles.Internal.secondOrderToString(zp.d2, significantDigits, name, normalized);

    // compute overall gain
    if normalized then
      gain :=zp.k*kn1*kn2/(kd1*kd2);
    else
      gain :=zp.k;
    end if;
    // Modelica.Utilities.Streams.print("gain = "+String(gain));

    // construct string for gain
    if gain <> 1.0 or gain == 1.0 and num_order == 0 then
      s := String(gain);
    end if;
    // Modelica.Utilities.Streams.print("s= "+s);

    // construct string for numerator
    if sn1 <> "" then
      if s == "" then
        s :=sn1;
      else
        s := s + "*" + sn1;
      end if;
    end if;
    if sn2 <> "" then
      if s == "" then
        s :=sn2;
      else
        s := s + "*" + sn2;
      end if;
    end if;

    // construct string for denominator
    if den_order <> 0 then
      s := s + " / ";
      if den_order > 1 then
        s := s + " ( ";
      end if;

      if sd1 <> "" then
        s := s + sd1;
      end if;

      if sd2 <> "" then
        if sd1 <> "" then
          s := s + "*";
        end if;
        s := s + sd2;
      end if;

      if den_order > 1 then
        s := s + " )";
      end if;
    end if;
  end 'String';

  encapsulated function p "Generate the transfer function p"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.ZerosAndPoles;

    output ZerosAndPoles zp(
      redeclare Real n1[1],
      redeclare Real n2[0,2],
      redeclare Real d1[0],
      redeclare Real d2[0,2]);
  algorithm
    zp.n1[1] := 0;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
p = ZerosAndPoles.<b>p</b>()
</pre></blockquote>

<h4>Description</h4>
<p>
Generate the complex Laplace variable p as a ZerosAndPoles transfer function. It can be used for generating like
</p>
<blockquote><pre>
ZerosAndPoles zp = p/(p^2 + p + 1)/(p + 1)
</pre></blockquote>
</html>"));
  end p;

  encapsulated package Analysis
    "Package of functions to analyse zeros-and-poles description represented by a ZerosAndPoles record"
    extends Modelica.Icons.Package;
    import Modelica;

    function analysis
      "Make a system analysis based on the poles and zeros of the system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Internal.AnalyseOptions;
      import Modelica_LinearSystems2.Internal.AnalyseOptions2;
      import Modelica_LinearSystems2.Internal.Eigenvalue;

      input ZerosAndPoles zp(uName="u", yName="y")
        "Transfer function of a system";

      input AnalyseOptions2 analyseOptions2 = Modelica_LinearSystems2.Internal.AnalyseOptions2(
               printControllability=false,
               printObservability=false);

      input String fileName="eigenvalues.html"
        "Name of html-file that contains eigenvalue table";

      input String systemName = ""
        "Name of system (used as heading in html file)";
      input String description = "" "Description of system (used in html file)";

    protected
      String dummyFileName = "dummy" + fileName;
      StateSpace ss=StateSpace(zp);

      AnalyseOptions analyseOptions=AnalyseOptions(
               plotEigenValues=analyseOptions2.plotEigenValues,
               plotInvariantZeros=analyseOptions2.plotInvariantZeros,
               plotStepResponse=analyseOptions2.plotStepResponse,
               plotFrequencyResponse=analyseOptions2.plotFrequencyResponse,
               printSystem = analyseOptions2.printSystem,
               printEigenValues=analyseOptions2.printEigenValues,
               printEigenValueProperties=analyseOptions2.printEigenValueProperties,
               printInvariantZeros=analyseOptions2.printInvariantZeros,
               printControllability=analyseOptions2.printControllability,
               printObservability=analyseOptions2.printObservability,
               headingEigenValues=analyseOptions2.headingEigenValues,
               headingInvariantzeros=analyseOptions2.headingInvariantzeros,
               headingStepResponse=analyseOptions2.headingStepResponse,
               headingFrequencyResponse=analyseOptions2.headingFrequencyResponse,
               dB_w=analyseOptions2.dB_w);

    algorithm
      assert(ZerosAndPoles.Analysis.denominatorDegree(zp) >= ZerosAndPoles.Analysis.numeratorDegree(zp),
        " Denominator polynominal of ZerosAndPoles object in function\"ZerosAndPoles.Analysis.analysis\"has to be of higher or equal order than numerator polynomial");
      Modelica.Utilities.Files.removeFile(fileName);
      Modelica.Utilities.Files.removeFile(dummyFileName);
      if analyseOptions.printSystem and size(ss.A,1) <= 50 then
        printSystem(
          zp,
          fileName,
          systemName,
          description);
        printSystem(
          zp,
          dummyFileName,
          systemName,
          description);
      end if;
      Modelica.Utilities.Streams.readFile(dummyFileName);
      analyseOptions.printSystem :=false;
      StateSpace.Analysis.analysis(
        ss=ss,
        analyseOptions=analyseOptions,
        fileName=fileName,
        systemName=systemName,
        description=description);

    public
     encapsulated function printSystem
        "Print the state space system in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.ZerosAndPoles;
        import Modelica_LinearSystems2;

       input ZerosAndPoles zp "transfer function to analyze";
       input String fileName="systemAnalysis.html"
          "File on which the zeros-and-poles transfer function is written in html format";
       input String systemName="ZerosAndPoles Transfer Function"
          "name of the system";
       input String description = ""
          "Description of system (used in html file)";
       input String format=".3g" "Format of numbers (e.g. \"20.8e\")";
      protected
       String st=String(zp);

     algorithm
       Modelica.Utilities.Files.removeFile(fileName);
       print("<html><body><br><br><p>\n<b>System report</b>\n</p>",fileName);
       print("<body><p><br> The system <b>" + systemName + "</b> is defined by</p>",fileName);
       print("G(p) = "+st, fileName);

       if description=="" then
         print("</table>", fileName);
       else
         print("</table>", fileName);
         print("<body><p><b>Description</b></p>",fileName);
         print(description, fileName);
       end if;

     end printSystem;

     annotation (__Dymola_interactive=true, Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
    end analysis;

   encapsulated function timeResponse
      "Calculate the time response of a zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

     extends Modelica_LinearSystems2.Internal.timeResponseMask2_zp;     // Input/Output declarations of time response functions
      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step;

     input Real x0[Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp)]=zeros(Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp))
        "Initial state vector";

    protected
     StateSpace ss=StateSpace(zp);

   algorithm
     (y,t,x_continuous) := StateSpace.Analysis.timeResponse(sc=ss, dt=dt, tSpan=tSpan, response=response, x0=x0);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = ZerosAndPoles.Analysis.<b>timeResponse</b>(zp, dt, tSpan, responseType, x0)
</pre></blockquote>

<h4>Description</h4>
<p>First, the ZerosAndPoles record is transformed into state space representation which is given to StateSpace.Analysis.timeResponse to calculate the time response of the state space system. The type of the time response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step &quot;Step response&quot;,
Ramp &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;
</pre></blockquote>
<p>
The state space system is transformed to a appropriate discrete state space system and, starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and x are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  p=Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=1/(p^2 + p + 1)

  Real Ts=0.1;
  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[2]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

algorithm
  (y,t,x):=Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse(zp,Ts,tSpan,response,x0);
//  y[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
</pre></blockquote>
</html>"));
   end timeResponse;

  encapsulated function impulseResponse "Calculate the impulse time response"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_zp;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse(
          zp=zp,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
          x0=zeros(Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp)));

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = ZerosAndPoles.Analysis.<b>impulseResponse</b>(zp, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>impulseResponse</b> calculates the time response of a ZerosAndPoles transfer function with impulse imput.
The system is first transformed zo a state space system, wich is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.impulseResponse(zp, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.timeResponse(zp, dt, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,ZerosAndPoles.Analysis.denominatorDegree(zp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.ZerosAndPoles zp=zp=1/(p^2 + p + 1)
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.ZerosAndPoles.Analysis.impulseResponse(zp,Ts,tSpan);
//  y[:,1,1]={0, 0.095, 0.18, 0.2553, 0.321}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.095, 0.18, 0.2553, 0.321}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse\">ZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end impulseResponse;

  encapsulated function stepResponse "Calculate the step time response"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_zp;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse(
          zp=zp,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
          x0=zeros(Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp)));

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = ZerosAndPoles.Analysis.<b>stepResponse</b>(zp, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>stepResponse</b> calculates the step response of a transfer function.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.stepResponse(zp, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.timeResponse(zp, dt, tSpan, response=Types.TimeResponse.Step, x0=fill(0,ZerosAndPoles.Analysis.denominatorDegree(zp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.ZerosAndPoles zp=zp=1/(p^2 + p + 1)
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.ZerosAndPoles.Analysis.stepResponse(zp,Ts,tSpan);
//  y[:,1,1]={0, 0.0048, 0.01867, 0.04, 0.0694}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0048, 0.01867, 0.04, 0.0694}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse\">ZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end stepResponse;

  encapsulated function rampResponse "Calculate the ramp time response"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_zp;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse(
          zp=zp,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
          x0=zeros(Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp)));

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = ZerosAndPoles.Analysis.<b>rampResponse</b>(zp, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rampResponse</b> calculates the time response of a transfer function for ramp imput u = t.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.rampResponse(zp, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.timeResponse(zp, dt, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,ZerosAndPoles.Analysis.denominatorDegree(zp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.ZerosAndPoles zp=zp=1/(p^2 + p + 1)
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.ZerosAndPoles.Analysis.rampResponse(zp,Ts,tSpan);
//  y[:,1,1]={0, 0.0002, 0.0012, 0.0042, 0.0096}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0002, 0.0012, 0.0042, 0.0096}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse\">ZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end rampResponse;

  encapsulated function initialResponse
      "Calculate the time response for given initial condition and zero inputs"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;

    input Real x0[:]=fill(0,0) "Initial state vector";

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_zp;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse(
          zp=zp,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
          x0=x0);

  annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = ZerosAndPoles.Analysis.<b>initialResponse</b>(zp, dt, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>initialResponse</b> calculates the time response of a state space system for given initial condition and zero inputs.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.initialResponse(x0,zp, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
ZerosAndPoles.Analysis.timeResponse(zp, dt, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.ZerosAndPoles zp=zp=1/(p^2 + p + 1)
  Real Ts=0.1;
  Real tSpan= 0.4;
  Real x0[2] = {1,1};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.ZerosAndPoles.Analysis.initialResponse(x0,zp,Ts,tSpan);
//  y[:,1,1]={1, 1.0903, 1.1616, 1.2151, 1.252}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 1.0903, 1.1616, 1.2151, 1.252}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.timeResponse\">ZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end initialResponse;

    encapsulated function numeratorDegree
      "Return numerator degree of a ZerosAndPoles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Integer result;
    algorithm
      result := size(zp.n1, 1) + 2*size(zp.n2, 1);
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = ZerosAndPoles.Analysis.<b>numeratorDegree</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>numeratorDegree</b> calculates the degree of the numerator polynomial constituted by the first and second order polynomials of the ZeroAndPoles numerator.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Real nDegree;

<b>algorithm</b>
  nDegree := ZerosAndPoles.Analysis.numeratorDegree(zp);
//  nDegree = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree\">ZerosAndPoles.Analysis.denominatorDegree</a>
</p>
</html>"));
    end numeratorDegree;

    encapsulated function denominatorDegree
      "Return denominator degree of a ZerosAndPoles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Integer result;
    algorithm
      result := size(zp.d1, 1) + 2*size(zp.d2, 1);
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = ZerosAndPoles.Analysis.<b>denominatorDegree</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial constituted by the first and second order polynomials of the ZeroAndPoles denominator.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Real dDegree;

<b>algorithm</b>
  dDegree := ZerosAndPoles.Analysis.denominatorDegree(zp);
//  dDegree = 2
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.numeratorDegree\">ZerosAndPoles.Analysis.numeratorDegree</a>
</p>
</html>"));
    end denominatorDegree;

    encapsulated function evaluate
      "Evaluate a ZerosAndPoles transfer function at a given value of p"
      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      input Complex p=Complex(0) "Complex value p where zp is evaluated";
      input Real den_min(min=0)=0 "|denominator(p)| is limited by den_min";
      output Complex y "= zp(p)";
    protected
      Complex num;
      Complex den;
      Real abs_den;
      Integer n1 = size(zp.n1,1);
      Integer n2 = size(zp.n2,1);
      Integer d1 = size(zp.d1,1);
      Integer d2 = size(zp.d2,1);
      Complex n[size(zp.n1,1)+size(zp.n2,1)];
      Complex d[size(zp.d1,1)+size(zp.d2,1)];
      Complex y2;
      Integer info;
    algorithm
      // Build numerator
      for i in 1:n1 loop
        n[i] :=ZerosAndPoles.Internal.'p+a'(p, zp.n1[i]);
      end for;
      for i in 1:n2 loop
        n[n1+i] := ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(p, zp.n2[i, :]);
      end for;

      // Build denominator
      for i in 1:d1 loop
        d[i] := ZerosAndPoles.Internal.'p+a'(p, zp.d1[i]);
      end for;
      for i in 1:d2 loop
        d[d1+i] :=ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(p, zp.d2[i, :]);
      end for;

      // Build value of transfer function
      (y2, info) :=Modelica_LinearSystems2.Internal.complexFraction(n, d);
      if info == 0 then
         y :=Complex(zp.k, 0)*y2;
      elseif info == 1 then
         y :=if zp.k >= 0 then y2 else -y2;
      else
         y :=y2;
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = ZerosAndPoles.Analysis.<b>evaluate</b>(zp, p, den_min=0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the ZerosAndPoles transfer function at a given (complex) value of p and returns the value G(p)=N(p)/D(p). The optional argument den_min with default 0 is used to guard against a division by zero.
</p>
<blockquote><pre>
<b>if</b> |(D(p))| >= den_min <b>then</b>
   G(p) = N(p) / D(p);
<b>elseif</b> D(p).re >= 0.0 <b>then</b>
   G(p) = N(p) / den_min
<b>else</b>
   G(p) = -N(p) / den_min
<b>end if</b>;
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Complex j = Modelica_LinearSystems2.Math.Complex.j();
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Complex result;

<b>algorithm</b>
  result := Modelica_LinearSystems2.ZerosAndPoles.Analysis.evaluate(zp, j+1);
//  result = 0.538462 - 0.307692j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>
</html>",     revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
    end evaluate;

    encapsulated function zerosAndPoles
      "Calculate zeros and poles of a ZerosAndPoles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Complex z[:]=fill(Complex(0, 0), size(zp.n1, 1) + 2*size(zp.n2, 1))
        "Zeros (Complex vector of numerator zeros)";
      output Complex p[:]=fill(Complex(0, 0), size(zp.d1, 1) + 2*size(zp.d2, 1))
        "Poles (Complex vector of denominator zeros)";
      output Real k "Constant multiplied with transfer function";

    protected
      Integer n_num1=size(zp.n1, 1);
      Integer n_num2=size(zp.n2, 1);
      Integer n_den1=size(zp.d1, 1);
      Integer n_den2=size(zp.d2, 1);
      Integer n_num=n_num1 + 2*n_num2;
      Integer n_den=n_den1 + 2*n_den2;
      Real re;
      Real im;
      Integer nz_real=ZerosAndPoles.Internal.numberOfRealZeros(zp.n1, zp.n2)
        "z[1:nz_real] are the real zeros";
      Integer np_real=ZerosAndPoles.Internal.numberOfRealZeros(zp.d1, zp.d2)
        "p[1:np_real] are the real poles";
      Real num_zeros1[nz_real];
      Real den_zeros1[np_real];
      Complex num_zeros2[:]=fill(Complex(0, 0), integer((n_num - nz_real)/2));
      Complex den_zeros2[:]=fill(Complex(0, 0), integer((n_den - np_real)/2));
      Integer n;
      Integer jj;

    algorithm
      (num_zeros1,num_zeros2) := ZerosAndPoles.Internal.roots(
          zp.n1,
          zp.n2,
          nz_real);
      (den_zeros1,den_zeros2) := ZerosAndPoles.Internal.roots(
          zp.d1,
          zp.d2,
          np_real);

      n := size(num_zeros1, 1);
      for i in 1:n loop
        z[i] := Complex(num_zeros1[i], 0);
      end for;

      jj := 1;
      for i in 1:size(num_zeros2, 1) loop
        z[n + jj] := num_zeros2[i];
        z[n + jj + 1] := num_zeros2[i].re-num_zeros2[i].im*j;
        jj := jj + 2;
      end for;

      n := size(den_zeros1, 1);
      for i in 1:n loop
        p[i] := den_zeros1[i]+0*j;
      end for;

      jj := 1;
      for i in 1:size(den_zeros2, 1) loop
        p[n + jj] := den_zeros2[i];
        p[n + jj + 1] := den_zeros2[i].re-den_zeros2[i].im*j;
        jj := jj + 2;
      end for;

      k := zp.k;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(z,p,k) = ZerosAndPoles.Analysis.<b>zerosAndPoles</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the zeros, poles and gain of a ZerosAndPoels transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

public
  output Complex z;
  output Complex p;
  output Real k;

<b>algorithm</b>
  (z,p,k)=Modelica_LinearSystems2.ZerosAndPoles.Analysis.zerosAndPoles(zp);
//  z = {-1}
//  p = {-0.5 + 0.866025j, -0.5 - 0.866025j}
//  k = 1
</pre></blockquote>
</html>"));
    end zerosAndPoles;

    function eigenValues
      "Calculate the eigen values of a linear zeros-and-poles transfer function and write them in a complex vector"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Complex eigval[:] "eigen values of the system";
    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      assert(ZerosAndPoles.Analysis.denominatorDegree(zp) >
        ZerosAndPoles.Analysis.numeratorDegree(zp),
        " Denominator polynominal of transfer function in function\"ZerosAndPoles.Analysis.eigenValues\"has to be of higher order than numerator polynomial");
      eigval := Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
eigenvalues = ZerosAndPoles.Analysis.<b>eigenValues</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvalues of the corresponding state space representation of a zeros-and-poles transfer function. The output is a complex vector containing the eigenvalues. Note, that the conversion of the transfer function does not result in a minimal state space system. Therefore also unobservable and uncontrollable eigenvalues will be calculated.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Complex eigenvalues[2];

<b>algorithm</b>
  eigenvalues = Modelica_LinearSystems2.ZerosAndPoles.Analysis.eigenValues(zp);
// eigenvalues = {-0.5 + j*sqrt(3)/2, -0.5 - j*sqrt(3)/2}
</pre></blockquote>
</html>"));
    end eigenValues;

    encapsulated function eigenVectors
      "Calculate the rigth eigenvectors of the corresponding state space system of a zeros-and-poles transfer function and write them columnwise in a matrix"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Complex;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      input Boolean onlyEigenvectors=true;
      output Real eigvec[:,:] "eigen values of the system";
      output Complex eigval[:] "eigen values of the system";
    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      assert(ZerosAndPoles.Analysis.denominatorDegree(zp) >
        ZerosAndPoles.Analysis.numeratorDegree(zp),
        " Denominator polynominal of transfer function in function\"ZerosAndPoles.Analysis.eigenVectors\"has to be of higher order than numerator polynomial");
      (eigvec,eigval) := StateSpace.Analysis.eigenVectors(ss=ss,
        onlyEigenvectors=onlyEigenvectors);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(eigenvectors, eigenvalues) = ZerosAndPoles.Analysis.<b>eigenVectors</b>(zp, onlyEigenvectors)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvectors and optionally (onlyEigenvectors=false) the eigenvalues of the corresponding state space system of a zeros-and-poles-transfer function. The output <tt>eigenvectors</tt> is a matrix with the same dimension as matrix <b>ss.A</b>. Just like in <a href=\"modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>, if the i-th eigenvalue has an imaginary part, then <tt>eigenvectors</tt>[:,i] is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the eigenvector of the i-th eigenvalue.<br>
The eigenvalues are returned as a complex vector <tt>eigenvalues</tt>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(2*p+2)/(p^2+2*p+2);

  Real eigenvectors[2,2];
  Complex eigenvalues[2];

<b>algorithm</b>
  (eigenvectors, eigenvalues) = Modelica_LinearSystems2.ZerosAndPoles.Analysis.eigenVectors(zp, true);
// eigenvectors = [(-0.4082), (-0.4082);
                    0.8165, 0]
// eigenvalues = {-1 + 1j, -1 - 1j}

          |-0.4082 -i0.4082 |         | -0.4082 + i0.4082 |
i.e. v1 = |                 |,   v2 = |                   |
          |     0.8165      |         |      0.8165       |
</pre></blockquote>
</html>"));
    end eigenVectors;

    encapsulated function invariantZeros
      "Compute invariant zeros of zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Complex;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";

      output Complex Zeros[:]= ZerosAndPoles.Analysis.zerosAndPoles(zp);

    algorithm

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zeros = ZerosAndPoles.Analysis.<b>invariantZeros</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes the invariant zeros of the corresponding state space representation of a zeros-and-poles transfer function. The output is a complex vector containing the eigenvalues. Note, that the conversion of the transfer function does not result in a minimal state space system. Therefore, also zeros equal to unobservable or uncontrollable eigenvalues will be computed.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Complex zeros[:];

<b>algorithm</b>
  zeros := Modelica_LinearSystems2.ZerosAndPoles.Analysis.invariantZeros(zp);
// zeros = {-1}

</pre></blockquote>
</html>",   revisions=
             "<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
    end invariantZeros;

    encapsulated function dcGain
      "Return steady state gain k (for a stable system: k = value of y at infinite time for a step input)"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica.Utilities.Streams.print;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Real k "Steady state gain";
      output Boolean finite = true
        "True, if k is finite, otherwise if k is infinite (k=Modelica.Constants.inf returned)";
    protected
      StateSpace ss=StateSpace(zp);
      Real K[1,1];
    algorithm
      (K, finite) := StateSpace.Analysis.dcGain(ss=ss);
      k :=K[1, 1];

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
          k = <b>dcGain</b>(zp);
(k, finite) = <b>dcGain</b>(zp);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the steady state gain <b>k</b> of a
ZerosAndPoles transfer function g(s), i.e. k = g(s=0).
For a stable transfer function, a step input u results
in the output y(t->t<sub>&infin;</sub>) = k.
</p>

<p>
If the transfer function has one or more zero poles, <b>k</b> is infinite.
In this case, the output argument <b>finite</b> = <b>false</b> and
<b>k</b> = Modelica.Constants.inf.
</p>
</html>"));
    end dcGain;

    encapsulated function isControllable
      "Check controllability of a zp-transfer-function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

      output Boolean controllable;
    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      controllable := StateSpace.Analysis.isControllable(ss=ss, method=method);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
controllable = ZerosAndPoles.Analysis.<b>isControllable</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function ZerosAndPoles.Analysis.<b>isControllable</b> checks the controllability
of a zeros-and-poles transfer function. Therefore, the transfer function is converted
into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isControllable\">StateSpace.Analysis.isControllable</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2 + 2*p +1);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean controllable;

<b>algorithm</b>
  controllable := Modelica_LinearSystems2.StateSpace.Analysis.isControllable(zp, method);
// controllable = true
</pre></blockquote>
</html>"));
    end isControllable;

    encapsulated function isObservable
      "Check observability of a zp-transfer-function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

        input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

        output Boolean observable;
    protected
        StateSpace ss=StateSpace(zp);

    algorithm
        observable := StateSpace.Analysis.isObservable(ss=ss, method=method);

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
observable = ZerosAndPoles.Analysis.<b>isObservable</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function ZerosAndPoles.Analysis.<b>isObservable</b> checks the observability of a zeros-and-poles transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isObservable\">StateSpace.Analysis.isObservable</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2 + 2*p +1);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean observable;

<b>algorithm</b>
  observable := Modelica_LinearSystems2.StateSpace.Analysis.isObservable(zp, method);
// observable = false
</pre></blockquote>
</html>"));
    end isObservable;

    encapsulated function isStabilizable
      "Check stabilizability of a zp-transfer-function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";

      output Boolean stabilizable;
    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      stabilizable := StateSpace.Analysis.isStabilizable(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
stabilizable = ZerosAndPoles.Analysis.<b>isStabilizable</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function ZerosAndPoles.Analysis.<b>isStabilizable</b> checks the Stabilizability of a zeros-and-poles transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isStabilizable\">StateSpace.Analysis.isStabilizable</a>.
The transfer function is stabilizable if all unstable poles are controllable.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p-1)/(p^2 - 2*p +1);

  Boolean stabilizable;

<b>algorithm</b>
   stabilizable := Modelica_LinearSystems2.ZerosAndPoles.Analysis.isStabilizable(zp);
// stabilizable = true
</pre></blockquote>
</html>"));
    end isStabilizable;

    encapsulated function isDetectable
      "Check detectability of a zp-transfer-function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";

      output Boolean detectable;

    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      detectable := StateSpace.Analysis.isDetectable(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
detectable = ZerosAndPoles.Analysis.<b>isDetectable</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function ZerosAndPoles.Analysis.<b>isDetectable</b> checks the Detectability of a zeros-and-poles transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isDetectable\">StateSpace.Analysis.isDetectable</a>. <br>
The transfer function is detectable if all unstable poles are observable.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p-1)/(p^2 - 2*p +1);

  Boolean detectable;

<b>algorithm</b>
  detectable := Modelica_LinearSystems2.ZerosAndPoles.Analysis.isDetectable(zp);
// detectable = false
</pre></blockquote>
</html>"));
    end isDetectable;

    encapsulated function controllabilityMatrix
      "Calculate the controllability matrix [B, A*B, ..., A^(n-1)*B] of a zp-transfer-function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Real om[:,:];

    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      om := StateSpace.Analysis.controllabilityMatrix(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = ZerosAndPoles.Analysis.<b>controllabilityMatrix</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the controllability matrix
</p>
<blockquote>
<b>Q</b> = [<b>B</b>, <b>A</b>*<b>B</b>, ..., <b>A</b>^(n-1)*<b>B</b>]
</blockquote>
<p>
of the system corresponding to state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
</pre></blockquote>
<p>
of a zeros and poles transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Real Q[2,2];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.ZerosAndPoles.Analysis.controllabilityMatrix(zp);
// Q = [0, 1, 1, -1]
</pre></blockquote>
</html>"));
    end controllabilityMatrix;

    encapsulated function observabilityMatrix
      "Calculate the observability matrix of zp-transfer-function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Real om[:,:];

    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      om := StateSpace.Analysis.observabilityMatrix(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = ZerosAndPoles.Analysis.<b>observabilityMatrix</b>(zp, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the observability matrix
</p>
<blockquote>
<b>Q</b> = [<b>C</b>; <b>C</b>*<b>A</b>; ...; <b>C</b>*<b>A</b>^(n-1)]
</blockquote>
<p>
of the system corresponding state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
</pre></blockquote>
<p>
of a zeros-and-poles transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Real Q[2,2];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.ZerosAndPoles.Analysis.observabilityMatrix(zp);
// Q = [1, 1, -1, 0]
</pre></blockquote>
</html>"));
    end observabilityMatrix;

  end Analysis;

  encapsulated package Design
    "Package of functions to design zeros-and-poles controllers and observers"
    extends Modelica.Icons.Package;
    import Modelica;

  encapsulated function filter
      "Generate a ZerosAndPoles transfer function from a filter description"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Utilities.Types;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input Modelica_LinearSystems2.Utilities.Types.AnalogFilter analogFilter=Modelica_LinearSystems2.Utilities.Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
      input Modelica_LinearSystems2.Utilities.Types.FilterType filterType=Modelica_LinearSystems2.Utilities.Types.FilterType.LowPass "Type of filter (LowPass/HighPass/BandPass)";
    input Integer order(min=1) = 2 "Order of filter";
    input Modelica.SIunits.Frequency f_cut=1/(2*Modelica.Constants.pi)
        "Cut-off frequency (default is w_cut = 1 rad/s)";
    input Real gain=1.0
        "Gain (= amplitude of frequency response at zero frequency)";
    input Real A_ripple(unit="dB") = 0.5
        "Pass band ripple (only for Chebyshev filter)";
    input Boolean normalized=true "True, if amplitude at f_cut = -3db*gain";
    input Modelica.SIunits.Frequency f_min=0
        "Band of normalized band pass/stop filter is f_min (-3db*gain) .. f_cut (-3db*gain)";

   /*
                               | size(n1,1)   | size(n2,1)       | size(d1,1)   | size(d2,1)
    ---------------------------+--------------+------------------+--------------+-------------------
    CriticalDamping - LowPass  |      0       |      0           |    order     |     0
    CriticalDamping - HighPass |    order     |      0           |    order     |     0
    CriticalDamping - BandPass |    order     |      0           |      0       |   order
    CriticalDamping - BandStop |      0       |    order         |      0       |   order
    ---------------------------+--------------+------------------+--------------+-------------------
    LowPass                    |      0       |      0           | mod(order,2) | integer(order/2)
    HighPass                   | mod(order,2) | integer(order/2) | mod(order,2) | integer(order/2)
    BandPass                   |     order    |      0           |      0       |   order
    BandStop                   |      0       |    order         |      0       |   order
 */
    output ZerosAndPoles filter(
      redeclare Real n1[if filterType == Types.FilterType.BandPass then order else
                        if filterType == Types.FilterType.HighPass then
                           (if analogFilter == Types.AnalogFilter.CriticalDamping then order else
                               mod(order, 2)) else 0],
      redeclare Real n2[if filterType == Types.FilterType.BandStop then order else
                        if filterType == Types.FilterType.HighPass and
                           analogFilter <> Types.AnalogFilter.CriticalDamping then integer(order/2) else 0,2],
      redeclare Real d1[if filterType == Types.FilterType.BandPass or
                           filterType == Types.FilterType.BandStop then 0 else
                        if analogFilter == Types.AnalogFilter.CriticalDamping then order else
                           mod(order, 2)],
      redeclare Real d2[if filterType == Types.FilterType.BandPass or
                           filterType == Types.FilterType.BandStop then order else
                        if analogFilter == Types.AnalogFilter.CriticalDamping then 0 else
                           integer(order/2),2]) "Filter transfer function";
    protected
    Integer n_num1=size(filter.n1, 1);
    Integer n_num2=size(filter.n2, 1);
    Integer n_den1=size(filter.d1, 1);
    Integer n_den2=size(filter.d2, 1);
    Integer n_num=n_num1 + 2*n_num2;
    Integer n_den=n_den1 + 2*n_den2;
    Real pi=Modelica.Constants.pi;
    Boolean evenOrder=mod(order, 2) == 0
        "True, if even filter order, otherwise uneven";

    Modelica.SIunits.Frequency f0 = if filterType == Types.FilterType.BandPass or
                                       filterType == Types.FilterType.BandStop then
                                                     sqrt(f_min*f_cut) else f_cut;
    Modelica.SIunits.AngularVelocity w_cut=2*pi*f0 "Cut-off angular frequency";
    constant Modelica.SIunits.AngularVelocity wOne = 1.0 "Just to make unit handling correct";
    Modelica.SIunits.AngularVelocity w_band = wOne*(f_cut - f_min) / f0;
    Real w_cut2 "= w_cut*w_cut";
    Real alpha=1.0 "Frequency correction factor";
    Real alpha2 "= alpha*alpha";

    Real alphax;

    Real epsilon "Ripple size";
    Real fac "arsinh(epsilon)";
    Real A2 "poleReal^2 + poleImag^2";
    Real A "Amplitude at w_cut";

    Real aux;
    Real k;
    Integer j;
    Real c;
    Real ww;
    ZerosAndPoles baseFilter;
  algorithm
    /* Compute filter coefficients of prototype low pass filter. If another filter
     characteristics is desired (e.g. high pass filter), it is derived
     from the low pass filter coefficients below
  */
    assert(f_cut > 0, "Cut-off frequency f_cut must be positive");
    baseFilter :=ZerosAndPoles.Internal.baseFilter(analogFilter,order,A_ripple,normalized);

    /* ==============================================================================
     Compute desired filter characteristics by transformation of low pass filter
  */
    if filterType == Types.FilterType.LowPass then
       filter :=baseFilter;

    elseif filterType == Types.FilterType.HighPass then
       /* The high pass filter is derived from the low pass filter by
        the transformation new(p) = 1/p
        1/(p + a)         -> 1/((1/p) + a) = (1/a)*p / (p + (1/a))
        1/(p^2 + a*p + b) -> 1/((1/p)^2 + a*(1/p) + b) = (1/b)*p^2 / (p^2 + (a/b)*p + 1/b)
     */
      assert(n_num1 == n_den1 and n_num2 == n_den2, "Internal error 1, should not occur");
      filter.k  := baseFilter.k;
      filter.n1 := zeros(n_num1);
      filter.n2 := zeros(n_num2, 2);

      for i in 1:n_den1 loop
        filter.d1[i] := 1/baseFilter.d1[i];
      end for;

      for i in 1:n_den2 loop
        filter.d2[i, 1] := baseFilter.d2[i, 1]/baseFilter.d2[i, 2];
        filter.d2[i, 2] := 1/baseFilter.d2[i, 2];
      end for;

    elseif filterType == Types.FilterType.BandPass then
      /* The band pass filter is derived from the low pass filter by
       the transformation new(p) = (p + 1/p)/w   (w = w_band = (f_max - f_min)/sqrt(f_max*f_min) )

       1/(p + a)         -> 1/(p/w + 1/p/w) + a)
                            = w*p / (p^2 + a*w*p + 1)

       1/(p^2 + a*p + b) -> 1/( (p+1/p)^2/w^2 + a*(p + 1/p)/w + b )
                            = 1 / ( p^2 + 1/p^2 + 2)/w^2 + (p + 1/p)*a/w + b )
                            = w^2*p^2 / (p^4 + 2*p^2 + 1 + (p^3 + p)a*w + b*w^2*p^2)
                            = w^2*p^2 / (p^4 + a*w*p^3 + (2+b*w^2)*p^2 + a*w*p + 1)

                            Assume the following description with PT2:
                            = w^2*p^2 /( (p^2 + p*(c/alpha) + 1/alpha^2)*
                                         (p^2 + p*(c*alpha) + alpha^2) )
                            = w^2*p^2 / ( p^4 + c*(alpha + 1/alpha)*p^3
                                              + (alpha^2 + 1/alpha^2 + c^2)*p^2
                                              + c*(alpha + 1/alpha)*p + 1 )

                            and therefore:
                              c*(alpha + 1/alpha) = a*w           -> c = a*w / (alpha + 1/alpha)
                                                                       = a*w*alpha/(1+alpha^2)
                              alpha^2 + 1/alpha^2 + c^2 = 2+b*w^2 -> equation to determine alpha
                              alpha^4 + 1 + a^2*w^2*alpha^4/(1+alpha^2)^2 = (2+b*w^2)*alpha^2
                              or z = alpha^2
                              z^2 + a^2*w^2*z^2/(1+z)^2 - (2+b*w^2)*z + 1 = 0
    */
      assert(n_num2 == 0 and n_den1 == 0, "Internal error 2, should not occur");
      assert(n_num1 == order and n_den2 == order, "Internal error 3, should not occur:");
      assert(f_min > 0 and f_min < f_cut, "Lower band pass frequency f_min must be > 0 and < f_cut");

      filter.n1 := zeros(n_num1);
      filter.k  := baseFilter.k;

      for i in 1:size(baseFilter.d1,1) loop
         filter.k := filter.k*w_band;
         filter.d2[i,1] := baseFilter.d1[i]*w_band;
         filter.d2[i,2] := 1.0;
      end for;

      for i in 1:size(baseFilter.d2,1) loop
         filter.k := filter.k*w_band*w_band;
         alpha := ZerosAndPoles.Internal.bandPassAlpha(baseFilter.d2[i,1], baseFilter.d2[i,2], w_band);
         c     := baseFilter.d2[i,1]*w_band / (alpha + 1/alpha);
         j     := size(baseFilter.d1,1) + 2*i - 1;
         filter.d2[j,   1] := c/alpha;
         filter.d2[j,   2] := 1/alpha^2;
         filter.d2[j+1, 1] := c*alpha;
         filter.d2[j+1, 2] := alpha^2;
      end for;

    elseif filterType == Types.FilterType.BandStop then
      /* The band stop filter is derived from the low pass filter by
       the transformation new(p) = w/( (p + 1/p) )   (w = w_band = (f_max - f_min)/sqrt(f_max*f_min) )

       1/(p + a)         -> 1/( w/(p + 1/p) ) + a)
                            = 1/a*(p^2 + 1) / (p^2 + (w/a)*p + 1)

       1/(p^2 + a*p + b) -> 1/( w^2/(p + 1/p)^2 + a*w/(p + 1/p) + b )
                            = 1/b*(p^2 + 1)^2 / (p^4 + a*w*p^3/b + (2+w^2/b)*p^2 + a*w*p/b + 1)

                            Assume the following description with PT2:
                            = 1/b*(p^2 + 1)^2 / ( (p^2 + p*(c/alpha) + 1/alpha^2)*
                                                  (p^2 + p*(c*alpha) + alpha^2) )
                            = 1/b*(p^2 + 1)^2 / (  p^4 + c*(alpha + 1/alpha)*p^3
                                                   + (alpha^2 + 1/alpha^2 + c^2)*p^2
                                                   + c*(alpha + 1/alpha)*p + 1 )

                            and therefore:
                              c*(alpha + 1/alpha) = a*w/b         -> c = a*w/(b*(alpha + 1/alpha))
                              alpha^2 + 1/alpha^2 + c^2 = 2+w^2/b -> equation to determine alpha
                              alpha^4 + 1 + (a*w/b*alpha^2)^2/(1+alpha^2)^2 = (2+w^2/b)*alpha^2
                              or z = alpha^2
                              z^2 + (a*w/b*z)^2/(1+z)^2 - (2+w^2/b)*z + 1 = 0

                            same as:  ww = w/b
                              z^2 + (a*ww*z)^2/(1+z)^2 - (2+b*ww)*z + 1 = 0  -> same equation as for BandPass

    */
      assert(n_num1 == 0 and n_den1 == 0, "Internal error 4, should not occur");
      assert(n_num2 == order and n_den2 == order, "Internal error 5, should not occur:");
      assert(f_min > 0 and f_min < f_cut, "Lower band pass frequency f_min must be > 0 and < f_cut");

      filter.k :=baseFilter.k;
      for i in 1:order loop
        filter.n2[i,1] := 0.0;
        filter.n2[i,2] := 1.0;
      end for;

      for i in 1:size(baseFilter.d1,1) loop
        filter.d2[i,1] := w_band/baseFilter.d1[i];
        filter.d2[i,2] := 1.0;
      end for;

      for i in 1:size(baseFilter.d2,1) loop
         ww    := w_band/baseFilter.d2[i,2];
         alpha := ZerosAndPoles.Internal.bandPassAlpha(baseFilter.d2[i,1], baseFilter.d2[i,2], ww);
         c     := baseFilter.d2[i,1]*ww / (alpha + 1/alpha);
         j     := size(baseFilter.d1,1) + 2*i - 1;
         filter.d2[j,   1] := c/alpha;
         filter.d2[j,   2] := 1/alpha^2;
         filter.d2[j+1, 1] := c*alpha;
         filter.d2[j+1, 2] := alpha^2;
      end for;

    else
      Streams.error("analogFilter (= " + String(analogFilter) + ") is not supported");
    end if;

    /* Transform filter to desired cut-off frequency ===================================

     Change filter coefficients according to transformation new(p) = p/w_cut
     Numerator  :     (p/w)^2 + a*(p/w) + b = (1/w^2)*(p^2 + (a*w)*p + b*w^2)
                                  (p/w) + a = (1/w)*(p + w*a)
     Denominator: 1/((p/w)^2 + a*(p/w) + b) = w^2/(p^2 + (a*w)*p + b*w^2)
                              1/((p/w) + a) = w/(p + w*a)
  */
    w_cut2 := w_cut*w_cut;
    filter.n1 := w_cut*filter.n1;
    filter.d1 := w_cut*filter.d1;
    filter.n2 := [w_cut*filter.n2[:, 1],w_cut2*filter.n2[:, 2]];
    filter.d2 := [w_cut*filter.d2[:, 1],w_cut2*filter.d2[:, 2]];

    /* Add gain ======================================================================= */
    if filterType == Types.FilterType.LowPass then
       /* A low pass filter does not have numerator polynomials and all coefficients
        of the denominator polynomial are guaranteed to be non-zero. It is then
        easy to compute the gain:
           1/(p + a)         -> a/(p + a)        , since g(0) = 1; k = a
           1/(p^2 + a*p + b) -> b/(p^2 + a*p + b), since g(0) = 1; k = b
     */
       k := 1.0;
       for i in 1:n_den1 loop
         k :=k*filter.d1[i];
       end for;
       for i in 1:n_den2 loop
         k :=k*filter.d2[i, 2];
       end for;
       filter.k := gain*k;

     elseif filterType == Types.FilterType.HighPass or
            filterType == Types.FilterType.BandStop then
       /* A high pass filter and a band stop filter have g(s->infinity) = 1
        and therefore filter.k = 1 is required, in ZerosAndPoles formulation
     */
       filter.k := gain;

    elseif filterType == Types.FilterType.BandPass then
       /* The gain due to the w-trasnformation must be added */
       filter.k := filter.k*gain*w_cut^(2*n_den2-n_num1);

    else
       Streams.error("analogFilter (= " + String(analogFilter) + ") is not supported");
    end if;

    annotation (Documentation(info="<html>

<h4>Syntax</h4>

<blockquote><pre>
zp = <b>filter</b>(analogFilter, filterType, order, f_cut, gain, A_ripple, normalized);
</pre></blockquote>

<h4>Description</h4>

<p>
This function constructs a ZerosAndPoles transfer function
description of low and high pass filters. For more details see also
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Literature\">[Tietze2002]</a>, pp. 815-852.
</p>
<p>
Typical frequency responses for the 4 supported low pass filter types
are shown in the next figure (this figure was generated with function
<a href=\"modelica://Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter2\">Examples.ZerosAndPoles.plotBodeFilter2</a>):
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4Filters.png\">
</p>
<p>
The step responses of the same low pass filters are shown in the next figure,
starting from a steady state initial filter with initial input = 0.2:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4FiltersStepResponse.png\">
</p>
<p>
Obviously, the frequency responses give a somewhat wrong impression
of the filter characteristics: Although Butterworth and Chebyshev
filters have a significantly steeper magnitude as the
CriticalDamping and Bessel filters, the step responses of
the latter ones are much better. This means for example, that
a CriticalDamping or a Bessel filter should be selected,
if a filter is mainly used to make a non-linear inverse model
realizable.
</p>

<p>
Typical frequency responses for the 4 supported high pass filter types
are shown in the next figure:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4Filters.png\">
</p>
<p>
The corresponding step responses of these high pass filters are
shown in the next figure:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4FiltersStepResponse.png\">
</p>
<p>
All filters are available in <b>normalized</b> (default) and non-normalized form.
In the normalized form, the amplitude of the filter transfer function
at the cutoff frequency is 3 dB. Note, when comparing the filters
of this function with other software systems, the setting of &quot;normalized&quot;
has to be selected appropriately. For example, the signal processing
toolbox of Matlab provides the filters in non-normalized form and
therefore a comparison makes only sense, if normalized = <b>false</b>
is set.


</p>

<h4>Example</h4>
<blockquote><pre>
   Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping;
   Integer order=2;
   Modelica.SIunits.Frequency f_cut=10;

   ZerosAndPoles zp_filter;

<b>algorithm</b>
    zp_filter=Modelica_LinearSystems2.ZerosAndPoles.Design.filter(
      order=order,
      f_cut=f_cut,
      analogFilter=analogFilter);

// zp_filter = 9530.93/( (p + 97.6265)^2 )
</pre></blockquote>
</html>"));
  end filter;

  end Design;

  encapsulated package Plot
    "Package of functions to plot zeros and poles description responses"
    extends Modelica.Icons.Package;
    import Modelica;

  encapsulated function polesAndZeros
    "Plot eigenvalues and or the zeros of a zeros-and-poles transfer function"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Plot;

    input ZerosAndPoles zp "Linear system in ZerosAndPoles form";
    input Boolean poles=true "= true, to plot the poles of zp" annotation(choices(checkBox=true));
    input Boolean zeros=true "= true, to plot the zeros of zp" annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
       defaultDiagram = Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

    protected
    Real eval[ZerosAndPoles.Analysis.denominatorDegree(zp),2];
    Real invZerosRe[ZerosAndPoles.Analysis.numeratorDegree(zp)];
    Real invZerosIm[ZerosAndPoles.Analysis.numeratorDegree(zp)];

    Complex invZeros[:];
    Complex poles2[:];

    Plot.Records.Curve curves[2];
    Integer i;
    Plot.Records.Diagram diagram2;
  algorithm
    (invZeros, poles2) := ZerosAndPoles.Analysis.zerosAndPoles(zp);

    for i in 1:size(invZeros, 1) loop
       invZerosRe[i] := invZeros[i].re;
       invZerosIm[i] := invZeros[i].im;
    end for;

    for i in 1:size(poles2, 1) loop
       eval[i,1] := poles2[i].re;
       eval[i,2] := poles2[i].im;
    end for;

    i :=0;
    if poles then
       i :=i + 1;
       curves[i] :=Plot.Records.Curve(
                          x=eval[:, 1],
                          y=eval[:, 2],
                          legend="poles",
                          autoLine=false,
                          linePattern=Plot.Types.LinePattern.None,
                          lineSymbol=Plot.Types.PointSymbol.Cross);
    end if;

    if zeros then
       i :=i + 1;
       curves[i] :=Plot.Records.Curve(
                          x=invZerosRe,
                          y=invZerosIm,
                          legend="zeros",
                          autoLine=false,
                          linePattern=Plot.Types.LinePattern.None,
                          lineSymbol=Plot.Types.PointSymbol.Circle);
    end if;

       diagram2 :=defaultDiagram;
       diagram2.curve :=curves[1:i];
       Plot.diagram(diagram2,device);

     annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>polesAndZeros</b>(zp);
   or
ZerosAndPoles.Plot.<b>polesAndZeros</b>(
  zp,
  poles=true,
  zeros=true,
  plot=true,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>());
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a pole-zero-map of the poles and zeros of a transfer function
in ZerosAndPoles format. The Boolean inputs
&quot;poles&quot; and &quot;zeros&quot; define what to plot. If Boolean input &quot;plot = true&quot;, the pole-zero-map
is plotted. If false, only the diagram is generated and returned as output argument.
The records &quot;defaultDiagram&quot; and &quot;device&quot; allow to set various layout options and the
size and location of the diagram on the screen.
</p>

<h4>Example</h4>
<p>
The example <a href=\"modelica://Modelica_LinearSystems2.Examples.ZerosAndPoles.plotPolesAndZeros\">
Modelica_LinearSystems2.Examples.ZerosAndPoles.plotPolesAndZeros</a>
defines a transfer functions as:
</p>

<blockquote><pre>
  TransferFunction s  = TransferFunction.s();
  TransferFunction tf = (s^3 + 4*s + 1)/(s^4 + 2*s^3 + 3*s^2 + 4*s);
  ZerosAndPoles    zp = ZerosAndPoles(tf);

  Modelica_LinearSystems2.ZerosAndPoles.Plot.polesAndZeros(zp=zp,
      defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
            heading=\"Poles and zeros of \" + String(tf)));
</pre></blockquote>

<p>
and results in
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/polesAndZerosZP.png\">
</blockquote>
</html>"));
  end polesAndZeros;

  encapsulated function bode
    "Plot ZerosAndPoles transfer function as bode plot"
    import Modelica;
    import Complex;
    import SI = Modelica.SIunits;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Internal;
    import Modelica_LinearSystems2.Utilities.Plot;

    input ZerosAndPoles zp "ZerosAndPoles transfer function to be plotted";
    input Integer nPoints(min=2) = 200 "Number of points";
    input Boolean autoRange=true
      "True, if abszissa range is automatically determined";
    input SI.Frequency f_min(min=0) = 0.1
      "Minimum frequency value, if autoRange = false"
      annotation(Dialog(enable=not autoRange));
    input SI.Frequency f_max(min=0) = 10
      "Maximum frequency value, if autoRange = false"
      annotation(Dialog(enable=not autoRange));

    input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
    input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

    extends Internal.PartialPlotFunction(
      defaultDiagram=Internal.DefaultDiagramBodePlot(
        heading="Bode plot: " + String(zp)));

    input Boolean Hz=true
      "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)"
      annotation(choices(checkBox=true));
    input Boolean dB=false
      "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))"
      annotation(choices(checkBox=true),Dialog(enable=magnitude));

    input Boolean onFile=false
      "= true, if frequency response is stored on file as matrix [f,A,phi]"
      annotation(choices(checkBox=true));
    input String fileName="frequencyResponse.mat"
      "If onFile=true, file on which the frequency response will be stored"
      annotation(Dialog(enable=onFile));
    input String matrixName=if Hz and not dB then "fHz_A_phiDeg" elseif
                               Hz and dB then "fHz_AdB_phiDeg" elseif
                               not Hz and dB then "f_AdB_phiDeg" else "f_A_phiDeg"
      "If onFile=true, Name of matrix on file"                                                                            annotation(Dialog(enable=onFile));
    protected
    SI.AngularVelocity w[nPoints];
    SI.Frequency f[nPoints];
    SI.Conversions.NonSIunits.Angle_deg phi[nPoints];
    Real A[nPoints];
    Real fAp[nPoints,if onFile then 3 else 0];
    Boolean OK;
    Complex c;
    Integer window=0;
    SI.Angle phi_old;
    Complex numZeros[:];
    Complex denZeros[:];
    Plot.Records.Curve curves[2];
    Integer i;
    Plot.Records.Diagram diagram2[2];
    Boolean success;
  algorithm
    // Determine frequency vector f
    if autoRange then
      (numZeros,denZeros) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
    else
      numZeros := fill(Complex(0), 0);
      denZeros := fill(Complex(0), 0);
    end if;
    f := Internal.frequencyVector(
      nPoints,
      autoRange,
      f_min,
      f_max,
      numZeros,
      denZeros,
      defaultDiagram.logX);

    // Compute magnitude/phase at the frequency points
    phi_old := 0.0;
    for i in 1:nPoints loop
      w[i] := SI.Conversions.from_Hz(f[i]);
      c := ZerosAndPoles.Analysis.evaluate(zp, Complex(0, w[i]), 1e-10);
      A[i] := Modelica.ComplexMath.'abs'(c);
      phi_old := Modelica.ComplexMath.arg(c, phi_old);
      phi[i] := SI.Conversions.to_deg(phi_old);

      // Convert to other units, if required
      if not Hz then
        f[i] := w[i];
      end if;
      if dB then
        A[i] := 20*log10(A[i]);
      end if;
    end for;

    // Plot computed frequency response
    diagram2 := fill(defaultDiagram, 2);
    i := 0;
    if magnitude then
      i := i + 1;
      curves[i] := Plot.Records.Curve(
        x=f,
        y=A,
        autoLine=true);
      diagram2[i].curve := {curves[i]};
      diagram2[i].yLabel := if dB then "magnitude [dB]" else "magnitude";
      if phase then
        diagram2[i].xLabel:="";
      end if;
      if dB then
        diagram2[i].logY := false;
      end if;
    end if;

    if phase then
      i := i + 1;
      curves[i] := Plot.Records.Curve(
        x=f,
        y=phi,
        autoLine=true);
      diagram2[i].curve := {curves[i]};
      diagram2[i].yLabel := "phase [deg]";
      diagram2[i].logY := false;
      if magnitude then
        diagram2[i].heading:="";
      end if;
    end if;

    if not Hz then
       diagram2[i].xLabel:="Angular frequency [rad/s]";
    end if;

    if magnitude and phase then
      Plot.diagramVector(diagram2, device);
    else
      Plot.diagram(diagram2[1], device);
    end if;

    if onFile then
      fAp :=[f,A,phi];
      Modelica.Utilities.Files.removeFile(fileName);
      success:=Modelica.Utilities.Streams.writeRealMatrix(fileName,matrixName,fAp,append=false);
      if success then
        Modelica.Utilities.Streams.print("... Frequency response stored on file \"" +
          Modelica.Utilities.Files.fullPathName(fileName) + "\"");
      end if;
    end if;

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>bode</b>(zp)
   or
ZerosAndPoles.Plot.<b>bode</b>(
  zp,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the bode-diagram of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(p^2 + 5*p + 7)/(p + 2)/(p + 3);

<b>algorithm</b>
  Modelica_LinearSystems2.ZerosAndPoles.Plot.bode(zp)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</blockquote>
</html>"));
  end bode;

  encapsulated function timeResponse
      "Plot the time response of a system represented by a transfer function. The response type is selectable"
      import Modelica;
      import Modelica_LinearSystems2;
  //    import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.ZerosAndPoles zp;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
    input Real x0[ZerosAndPoles.Analysis.denominatorDegree(zp)]=zeros(
        ZerosAndPoles.Analysis.denominatorDegree(zp)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="time response of  zp = "
           + String(zp)));

    protected
    Plot.Records.Curve curve;
    Plot.Records.Diagram diagram2;
    Real y[:,1,1] "Output response";
    Real t[:] "Time vector: (number of samples)";

  algorithm
    (y,t) := ZerosAndPoles.Analysis.timeResponse(
      zp,
      dt,
      tSpan,
      response,
      x0);

    curve := Plot.Records.Curve(
      x=t,
      y=y[:, 1, 1],
      legend="y",
      autoLine=true);
    diagram2 := defaultDiagram;
    diagram2.curve := {curve};

    Plot.diagram(diagram2, device);
    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>timeResponse</b>(zp);
   or
ZerosAndPoles.Plot.<b>timeResponse</b>(
  zp,
  dt,
  tSpan,
  response,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the time response of a transfer function. The character of the time response if defined by the input
<a href=\"modelica://Modelica_LinearSystems2.Types.TimeResponse\">response</a>, i.e. Impulse, Step, Ramp, or Initial.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(p + 1)/(p^2 + 5*p + 12);

  Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
  Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse(zp, dt=0.02, tSpan=3, response=response)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/timeResponseZP.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end timeResponse;

  encapsulated function impulse "Impulse response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.ZerosAndPoles zp
        "zeros-and-poles transfer function";
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

    input Real x0[ZerosAndPoles.Analysis.denominatorDegree(zp)]=zeros(ZerosAndPoles.Analysis.denominatorDegree(zp))
        "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response of  zp = "
           + String(zp)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "type of time response";
  algorithm
    Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse(
      zp=zp,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>impulse</b>(zp)
   or
ZerosAndPoles.Plot.<b>impulse</b>(
  zp,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the impulse response of a zeros-and-poles transfer function. It is based on
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(p + 1)/(p^2 + 5*p + 12);

<b>algorithm</b>
   Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse(zp, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/impulseResponseZP.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end impulse;

  encapsulated function step "Step response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.ZerosAndPoles zp;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
    input Real x0[ZerosAndPoles.Analysis.denominatorDegree(zp)]=zeros(
        ZerosAndPoles.Analysis.denominatorDegree(zp)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  zp = "
           + String(zp)));

  algorithm
   Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse(
      zp=zp,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

  equation

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>step</b>(zp)
   or
ZerosAndPoles.Plot.<b>step</b>(
  zp,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the step response of a zeros-and-poles transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(p + 1)/(p^2 + 5*p + 12);

<b>algorithm</b>
  Modelica_LinearSystems2.ZerosAndPoles.Plot.step(zp, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/stepResponseZP.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end step;

  encapsulated function ramp "Ramp response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.ZerosAndPoles zp;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";
    input Real x0[ZerosAndPoles.Analysis.denominatorDegree(zp)]=zeros(
        ZerosAndPoles.Analysis.denominatorDegree(zp)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  zp = "
           + String(zp)));

  algorithm
   Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse(
      zp=zp,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>ramp</b>(zp)
   or
ZerosAndPoles.Plot.<b>ramp</b>(
  zp,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the ramp response of a zeros-and-poles transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp =(2*p^2 + 7*p + 13)/(p + 1)/(p^2 + 5*p + 12);

<b>algorithm</b>
  Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp(zp)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/rampResponseZP.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end ramp;

  encapsulated function initialResponse "Initial condition response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.ZerosAndPoles zp;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";
    input Real y0 "Initial output (for initial condition plot)";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Initial response of  zp = "
           + String(zp) + "  with y0 = " + String(y0)));

    protected
    Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(zp);
    Real x0[ZerosAndPoles.Analysis.denominatorDegree(zp)]=
        Modelica.Math.Matrices.equalityLeastSquares(
        ss.A,
        fill(0, size(ss.B, 1)),
        ss.C,
        vector(y0)) "Initial state vector (for initial condition plot)";
  algorithm

    Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse(
          zp=zp,
          dt=dt,
          tSpan=tSpan,
          response=response,
          x0=x0,
          defaultDiagram=defaultDiagram,
          device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ZerosAndPoles.Plot.<b>initialResponse</b>(zp)
   or
ZerosAndPoles.Plot.<b>initialResponse</b>(
  zp,
  dt,
  tSpan,
  y0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the initial response, i.e. the zeros input response of a zeros and poles transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp = (p + 1)/(p^2 + 5*p + 12);
  Real y0=1;

<b>algorithm</b>
   Modelica_LinearSystems2.ZerosAndPoles.Plot.initialResponseZP(zp, y0=y0, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/ZerosAndPoles/initialResponseZP.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp\">ramp</a>
</p>
</html>"));
  end initialResponse;

  end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of ZerosAndPoles data record"
    extends Modelica.Icons.Package;
    import Modelica;

    function toTransferFunction
      "Generate a TransferFunction data record from a ZerosAndPoles data record"
      //encapsulated function fromZerosAndPoles
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.ZerosAndPoles;
      //import Modelica_LinearSystems2.Internal;
      import Complex;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output TransferFunction tf(
        redeclare Real n[2*size(zp.n2,1)+size(zp.n1,1)+1],
        redeclare Real d[2*size(zp.d2,1)+size(zp.d1,1)+1]);

    protected
      Real k;
      Complex z[:];
      Complex p[:];
      Polynomial pn;
      Polynomial pd;
    algorithm
      (z,p,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
      pn := Polynomial(z)*Polynomial(k);
      pd := Polynomial(p);
      tf.n := pn.c;
      tf.d := pd.c;
      tf.uName := zp.uName;
      tf.yName := zp.yName;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tf = ZerosAndPoles.Conversion.toStateSpace<b>toTransferFunction</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a TransferFunction record
</p>
<blockquote><pre>
      n(s)     b0 + b1*s + ... + bn*s^n
tf = ------ = --------------------------
      d(s)     a0 + a1*s + ... + an*s^n
</pre></blockquote>
<p>
from a ZerosAndPoles record representated by first and second order numerator and denominator polynomials. The poles and zeros and the gain <tt>k</tt> are computed (<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs in the TransferFunction constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp = 1/(p + 3)/(p + 1)

<b>algorithm</b>
  tf:=Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction(zp);
//  tf = 1/( s^2 + 4*s + 3 )
</pre></blockquote>
</html>"));
    end toTransferFunction;

    encapsulated function toTransferFunctionMIMO
      "Generate a transfer function matrix  from zeros-and-poles transfer function matrix"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.TransferFunction;

      input ZerosAndPoles zp[:,:] "ZerosAndPoles transfer function of a system";

      output TransferFunction tf[size(zp, 1),size(zp, 2)];

    protected
      Integer ny=size(zp, 1);
      Integer nu=size(zp, 2);

    algorithm
      for iy in 1:ny loop
        for iu in 1:nu loop
          tf[iy, iu] := ZerosAndPoles.Conversion.toTransferFunction(zp[iy, iu]);
        end for;
      end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tf = ZerosAndPoles.Conversion.toStateSpace<b>toTransferFunctionMIMO</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Converts a matrix of ZerosAndPoles transfer functions denoted by the product of first and second order numerator and denominator polynomials into a matrix of transfer functions represented by (usual) numerator and denominator polynomial. The function repetitively uses <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction\">toTransferFunction</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp = [1/(p + 2)/(p + 1);p/(p + 1)/(p + 1)]

<b>algorithm</b>
  tf := Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction(zp);
//  tf = [1/( (p + 1)*(p + 2) ); p/( (p + 1)^2 )]
</pre></blockquote>
</html>"));
    end toTransferFunctionMIMO;

    function toMatrices
      "Convert a ZerosAndPoles object into the matrices A, B, C of a StateSpace"
     //encapsulated function fromZerosAndPoles
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output Real ABCD[ZerosAndPoles.Analysis.denominatorDegree(zp)+1,
          ZerosAndPoles.Analysis.denominatorDegree(zp)+1];

    protected
      Real ssA[ZerosAndPoles.Analysis.denominatorDegree(zp),ZerosAndPoles.Analysis.denominatorDegree(zp)]
        "system matrix of partial 2nd order system";
      Real ssB[ZerosAndPoles.Analysis.denominatorDegree(zp),1]
        "input matrix of partial 2nd order system";
      Real ssC[1,ZerosAndPoles.Analysis.denominatorDegree(zp)]
        "output matrix of partial 2nd order system";
      Real ssD[1,1] "feedthrough matrix of partial 2nd order system";

      Real A[2,2] "system matrix of partial 2nd order system";
      Real B[2,1] "input matrix of partial 2nd order system";
      Real C[1,2] "output matrix of partial 2nd order system";
      Real D[1,1] "feedthrough matrix of partial 2nd order system";
      Real a "system 'matrix' of partial 1st order system";
      Real b "input 'matrix' of partial 1st order system";
      Real c "output 'matrix' of partial 1st order system";
      Real d "feedthrough 'matrix' of partial 1st order system";
      Real eps=1e-6;
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
      Integer ili=if max(n_den2, n_num2) > 0 then i_d else  max(2, i_d);

      Real num[nx,2]=[zp.n2; [zp.n1,zeros(n_num1)]; zeros(max(0,nx - n_num2 - n_num1), 2)]
        "Numerator matrix, in order that indices are defined in all situations in all if clauses";
      Real den[nx,2]=[zp.d2; [zp.d1,zeros(n_den1)]; zeros(max(0,nx - n_den2 - n_den1), 2)]
        "Denominator matrix, in order that indices are defined in all situations in all if clauses";
      Real k[i_k + n_den1](each fixed=false)
        "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
      Real k_total;

      Boolean dZero=true;

     //ZerosAndPoles zp2;

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
              k[i] := StateSpace.Internal.scaleFactor2(
                  num[i, 1],
                  num[i, 2],
                  den[i, 1],
                  den[i, 2],eps);
            elseif 2*(i - n_num2) <= n_num1 then
                // State space system in form (1) with 2 first order numerator polynomials
              k[i] := StateSpace.Internal.scaleFactor2(
                  num[max(1,2*(i - n_num2)-1), 1] + num[max(1,2*(i - n_num2)), 1],
                  num[max(1,2*(i - n_num2)-1), 1]*num[max(1,2*(i - n_num2)), 1],
                  den[i, 1],
                  den[i, 2],eps);
            elseif  2*(i-n_num2) -1== n_num1 then
                // State space system in form (2) with 1 first order numerator polynomial
              k[i] := StateSpace.Internal.scaleFactor2(
                  1,
                  num[2*i-n_num2-1, 1],
                  den[i, 1],
                  den[i, 2],eps);
            else
                // State space system in form (3)
              k[i] := StateSpace.Internal.scaleFactor2(
                  1,
                  1,
                  den[i, 1],
                  den[i, 2],eps);
            end if;
          else
             // State space system in form (1) with 2 first order denominator polynomials
            k[i] := StateSpace.Internal.scaleFactor2(
                num[i, 1],
                num[i, 2],
                den[max(1,2*(i - n_den2)-1), 1] + den[max(1,2*(i - n_den2)), 1],
                den[max(1,2*(i - n_den2)-1), 1]*den[max(1,2*(i - n_den2)), 1],eps);
          end if;
        end for;

        for i in i_d:n_den1 loop
          // State space systems of order 1
          if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
             // State space system in form (4)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
              n_num2) + i), 1], den[n_den2 + i, 1],eps);
          elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
             // State space system in form (4)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
              1], den[n_den2 + i, 1],eps);
          else
             // State space system in form (5)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(1, den[n_den2 + i, 1],eps);
          end if;
        end for;

        k_total := zp.k/product(k);

        ssA := zeros(nx, nx);
        ssB := zeros(nx, 1);
        ssC := zeros(1, nx);
        ssD := zeros(1, 1);

     // Calculation of matrices A, B, C, D
     //first elements of A, B, C and D
        if max(n_den2, n_num2) > 0 then
          A[1, :] := {0,1};
          B[1, 1] := 0;
          // Construct state space systems of order 2
          if 1 <= n_den2 then
            A[2, :] := {-den[1, 2],-den[1, 1]};
            B[2, 1] := if abs(den[1, 2])>eps and abs(num[1, 2])>eps then den[1, 2] else 1;

            if 1 <= n_num2 then
              // State space system in form (1)
              C := if abs(den[1, 2])>eps and abs(num[1, 2])>eps then k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
              D := [k[1]];
              dZero := false;
           elseif 1 - n_num2 + 1 <= n_num1 then
            // State space system in form (1) with 2 first order numerator polynomials
              B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then   k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
              D := [k[1]];
              dZero := false;
           elseif 1 - n_num2 == n_num1 then
                   // State space system in form (2) with 1 first order numerator polynomial
              B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then k[1]*[num[1, 1],1]/den[1, 2] else k[1]*[num[1, 1],1];
              D := [0];
              dZero := dZero and true;
            else
                   // State space system in form (3)
              B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then  k[1]*[1,0]/den[1, 2] else k[1]*[1,0];
              D := [0];
              dZero := dZero and true;
            end if;
          else
            // State space system in form (1) with 2 first order denominator polynomials
            A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
            B[2, 1] := if abs(den[1, 1]*den[2, 1])>eps then den[1, 1]*den[2, 1] else 1;
            C := if abs(den[1, 1]*den[2, 1])>eps then  k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/den[1, 1]/den[2, 1] else k[1]*[num[1, 2],num[1, 1] - den[2, 1]];
            D := [k[1]];
            dZero := false;
          end if;
          ssA[1:2, 1:2] := A;
          ssB[1:2, 1] := vector(B);
          ssC[1, 1:2] := vector(C);
          ssD := D;

        else

       // Construct state space systems of order 1
          a := -den[1, 1];

          if 1 <= n_num1 then
                // State space system in form (4)
            b := if abs(den[1, 1])>eps then den[1,1] else num[1,1];
            c := if abs(den[1, 1])>eps then k[1]*(num[1, 1] - den[1, 1])/den[1, 1] else k[1];
            d := k[1];
            dZero := false;
          else
                // State space system in form (5)
            b := if abs(den[1, 1])>eps then den[1,1] else if n_num1>0 then num[1,1] else 1;
            c := if abs(den[1, 1])>eps then k[1]/den[1, 1] else k[1];
            d := 0;
            dZero := dZero and true;
          end if;
          ssA[1, 1] := a;
          ssB[1, 1] := b;
          ssC[1, 1] := c;
          ssD[1, 1] := d;

        end if;
     /// for i=2 to degree(system)
        A[1, :] := {0,1};
        B[1, 1] := 0;
        for i in 2:max(n_den2, n_num2) loop
             // Construct state space systems of order 2
          if i <= n_den2 then
            A[2, :] := {-den[i, 2],-den[i, 1]};
            B[2, 1] := if abs(den[i, 2])>eps and abs(num[i, 2])>eps then den[i, 2] else 1;
            if i <= n_num2 then
                   // State space system in form (1)
              C := if abs(den[i, 2])>eps and abs(num[i, 2])>eps then k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
              D := [k[i]];
              dZero := false;
            elseif 2*(i - n_num2) <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
              C := if abs(den[i, 2])>eps and abs(num[2*i-n_num2-1, 2])>eps then k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]];
              D := [k[i]];
              dZero := false;
            elseif 2*(i-n_num2) -1== n_num1 then
            // State space system in form (2) with 1 first order numerator polynomial
            B[2, 1] := if abs(den[i, 2])>eps then den[i, 2] else 1;
              C := if abs(den[i, 2])>eps then k[i]*[num[2*i-n_num2-1, 1],1]/den[i, 2] else k[i]*[num[2*i-n_num2-1, 1],1];
              D := [0];
              dZero := dZero and true;
            else
              // State space system in form (3)
              B[2, 1] := if abs(den[i, 2])>eps then den[i, 2] else 1;
              C := if abs(den[i, 2])>eps then k[i]*[1,0]/den[i, 2] else k[i]*[1,0];
              D := [0];
              dZero := dZero and true;
            end if;

          else
                // State space system in form (1) with 2 first order denominator polynomials
            A[2, :] := {-(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),-(den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])};
            B[2, 1] := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1] else 1;
            C := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/den[max(2*(i-n_den2)-1,i), 1]/den[max(2*(i-n_den2),i), 1] else k[i]*[num[max(2*(i-n_num2),i), 2],num[max(2*(i-n_num2),i), 1] - den[max(2*(i-n_den2),i), 1]];
            D := [k[i]];
            dZero := false;
          end if;
          ssA[2*i, 1:2*i - 2] := B[2, 1]*ssC[1, 1:2*i - 2];
          ssA[2*i - 1:2*i, 2*i - 1:2*i] := A;
          ssB[2*i, 1] := if dZero then 0 else B[2, 1]*ssD[1, 1];
          ssC[1, 1:2*i - 2] := if dZero then fill(0, 2*i - 2) else D[1, 1]*ssC[
            1, 1:2*i - 2];
          ssC[1, 2*i - 1:2*i] := vector(C);
          ssD := D*ssD;
        end for;
     //  for i in max(2,i_d):n_den1 loop
        for i in ili:n_den1 loop
             // Construct state space systems of order 1
          a := if abs(den[n_den2 + i, 1])>eps then -den[n_den2 + i, 1] else 0.0;

          if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
                // State space system in form (4)

            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]*(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/den[n_den2 + i, 1] else 1.0;
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1];
            d := k[i_k + i];
            dZero := false;
          elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
                // State space system in form (4)
            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]*(num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/den[n_den2 + i, 1] else 1.0;
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else num[max(1, n_num2 + i - i_d + 1), 1];
            d := k[i_k + i];
            dZero := false;
          else
                // State space system in form (5)
            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]/den[n_den2 + i, 1] else k[i_k + i];
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else 1;
            d := 0;
            dZero := dZero and true;
          end if;
          ssA[2*n_den2 + i, 1:2*n_den2 + i - 1] := b*ssC[1, 1:2*n_den2 + i - 1];
          ssA[2*n_den2 + i, 2*n_den2 + i] := a;
          ssB[2*n_den2 + i, 1] := if dZero then 0 else b*ssD[1, 1];
          ssC[1, 1:2*n_den2 + i - 1] := if dZero then fill(0, 2*n_den2 + i - 1) else
                  d*ssC[1, 1:2*n_den2 + i - 1];
          ssC[1, 2*n_den2 + i] := c;
          ssD := if dZero then [0] else d*ssD;

        end for;

        ssC := k_total*ssC;
        ssD := k_total*ssD;
      else
        ABCD := [zp.k];
      end if;

      ABCD := [ssA,ssB; ssC,ssD];

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ABCD = ZerosAndPoles.Conversion.toStateSpace<b>toStateSpace</b>(zp)
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
G(p) = ------- = --- * ----------- = k * F(p),   k = a/b
        b + p     b       b + p
</pre></blockquote>
<p>
for first order systems respectively.
</p>
<p>
The complete system is now considered as the series connections of all
the single unit gain transfer functions and an overall gain k with
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
          b2 + b1*p + p^2                 y = c*((a2-b2)*x1 + (a1-b1)*x2 + u),  c = b2/a2
&nbsp;
(2)
             p + a                 der(x1) = x2
    y = ---------------- * u  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                  y = k*(a1/b2*x1 +x2/b2),  c = b2/a
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
  ABCD := Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace(zp);
// ssA = [0, 1; -1, -1],
// ssB = [0; 1],
// ssC = [1, 1],
// ssD = [0],
</pre></blockquote>
</html>"));
    end toMatrices;

    function toStateSpace
      "Transform a ZerosAndPoles object into a StateSpace object"
     //encapsulated function fromZerosAndPoles
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

      input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
      output StateSpace ss(
        redeclare Real A[ZerosAndPoles.Analysis.denominatorDegree(zp),
          ZerosAndPoles.Analysis.denominatorDegree(zp)],
        redeclare Real B[ZerosAndPoles.Analysis.denominatorDegree(zp),1],
        redeclare Real C[1,ZerosAndPoles.Analysis.denominatorDegree(zp)],
        redeclare Real D[1,1]) "Transfer function in StateSpace SISO form";

    protected
      Real A[2,2] "System matrix of partial 2nd order system";
      Real B[2,1] "Input matrix of partial 2nd order system";
      Real C[1,2] "Output matrix of partial 2nd order system";
      Real D[1,1] "Feedthrough matrix of partial 2nd order system";
      Real a "System 'matrix' of partial 1st order system";
      Real b "Input 'matrix' of partial 1st order system";
      Real c "Output 'matrix' of partial 1st order system";
      Real d "Feedthrough 'matrix' of partial 1st order system";
      Real eps = 1e-6;
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
      Integer ili=if max(n_den2, n_num2) > 0 then i_d else  max(2, i_d);

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
              k[i] := StateSpace.Internal.scaleFactor2(
                  num[i, 1],
                  num[i, 2],
                  den[i, 1],
                  den[i, 2],eps);
            elseif 2*(i - n_num2) <= n_num1 then
                // State space system in form (1) with 2 first order numerator polynomials
              k[i] := StateSpace.Internal.scaleFactor2(
                  num[2*(i - n_num2)-1, 1] + num[2*(i - n_num2), 1],
                  num[2*(i - n_num2)-1, 1]*num[2*(i - n_num2), 1],
                  den[i, 1],
                  den[i, 2],eps);
            elseif  2*(i-n_num2) -1== n_num1 then
                // State space system in form (2) with 1 first order numerator polynomial
              k[i] := StateSpace.Internal.scaleFactor2(
                  1,
                  num[2*i-n_num2-1, 1],
                  den[i, 1],
                  den[i, 2],eps);
             else
                // State space system in form (3)
              k[i] := StateSpace.Internal.scaleFactor2(
                  1,
                  1,
                  den[i, 1],
                  den[i, 2],eps);
            end if;
          else
             // State space system in form (1) with 2 first order denominator polynomials
            k[i] := StateSpace.Internal.scaleFactor2(
                num[i, 1],
                num[i, 2],
                den[2*(i - n_den2)-1, 1] + den[2*(i - n_den2), 1],
                den[2*(i - n_den2)-1, 1]*den[2*(i - n_den2), 1],eps);
          end if;
        end for;

        for i in i_d:n_den1 loop
          // State space systems of order 1
          if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
             // State space system in form (4)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
              n_num2) + i), 1], den[n_den2 + i, 1],eps);
          elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
             // State space system in form (4)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
              1], den[n_den2 + i, 1],eps);
          else
             // State space system in form (5)
            k[i_k + i] := StateSpace.Internal.scaleFactor1(1, den[n_den2 + i, 1],eps);
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
          A[1, :] := {0,1};
          B[1, 1] := 0;
          // Construct state space systems of order 2
          if 1 <= n_den2 then
            A[2, :] := {-den[1, 2],-den[1, 1]};
            B[2, 1] := if abs(den[1, 2])>eps and abs(num[1, 2])>eps then den[1, 2] else 1;

            if 1 <= n_num2 then
             // State space system in form (1)
              C := if abs(den[1, 2])>eps and abs(num[1, 2])>eps then k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
              D := [k[1]];
              dZero := false;
           elseif 1 - n_num2 + 1 <= n_num1 then
            // State space system in form (1) with 2 first order numerator polynomials
            B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then   k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/den[1, 2] else k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
              D := [k[1]];
              dZero := false;
           elseif 1 - n_num2 == n_num1 then
            // State space system in form (2) with 1 first order numerator polynomial
              B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then k[1]*[num[1, 1],1]/den[1, 2] else k[1]*[num[1, 1],1];
              D := [0];
              dZero := dZero and true;
            else
                   // State space system in form (3)
              B[2, 1] := if abs(den[1, 2])>eps then den[1, 2] else 1;
              C := if abs(den[1, 2])>eps then  k[1]*[1,0]/den[1, 2] else k[1]*[1,0];
              D := [0];
              dZero := dZero and true;
            end if;
          else
            // State space system in form (1) with 2 first order denominator polynomials
            A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
            B[2, 1] := if abs(den[1, 1]*den[2, 1])>eps then den[1, 1]*den[2, 1] else 1;
            C := if abs(den[1, 1]*den[2, 1])>eps then  k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/den[1, 1]/den[2, 1] else k[1]*[num[1, 2],num[1, 1] - den[2, 1]];
            D := [k[1]];
            dZero := false;
          end if;
          ss.A[1:2, 1:2] := A;
          ss.B[1:2, 1] := vector(B);
          ss.C[1, 1:2] := vector(C);
          ss.D := D;

        else

       // Construct state space systems of order 1
          a := -den[1, 1];

          if 1 <= n_num1 then
            // State space system in form (4)
            b := if abs(den[1, 1])>eps then den[1,1] else num[1,1];
            c := if abs(den[1, 1])>eps then k[1]*(num[1, 1] - den[1, 1])/den[1, 1] else k[1];
            d := k[1];
            dZero := false;
          else
           // State space system in form (5)
            b := if abs(den[1, 1])>eps then den[1,1] else if n_num1>0 then num[1,1] else 1;
            c := if abs(den[1, 1])>eps then k[1]/den[1, 1] else k[1];
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
            B[2, 1] := if abs(den[i, 2])>eps and abs(num[i, 2])>eps then den[i, 2] else 1;

            if i <= n_num2 then
                   // State space system in form (1)

              C := if abs(den[i, 2])>eps and abs(num[i, 2])>eps then k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
              D := [k[i]];
              dZero := false;

            elseif 2*(i - n_num2) <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
              C := if abs(den[i, 2])>eps and abs(num[2*i-n_num2-1, 2])>eps then k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]]/den[i, 2] else k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]];
              D := [k[i]];
              dZero := false;

            elseif 2*(i-n_num2) -1== n_num1 then
            // State space system in form (2) with 1 first order numerator polynomial
            B[2, 1] := if abs(den[i, 2])>eps then den[i, 2] else 1;
              C := if abs(den[i, 2])>eps then k[i]*[num[2*i-n_num2-1, 1],1]/den[i, 2] else k[i]*[num[2*i-n_num2-1, 1],1];
              D := [0];
              dZero := dZero and true;
            else
              // State space system in form (3)
              B[2, 1] := if abs(den[i, 2])>eps then den[i, 2] else 1;
              C := if abs(den[i, 2])>eps then k[i]*[1,0]/den[i, 2] else k[i]*[1,0];
              D := [0];
              dZero := dZero and true;
            end if;

          else
                // State space system in form (1) with 2 first order denominator polynomials
            A[2, :] := {-(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),-(den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])};
            B[2, 1] := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1] else 1;
            C := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/den[max(2*(i-n_den2)-1,i), 1]/den[max(2*(i-n_den2),i), 1] else k[i]*[num[max(2*(i-n_num2),i), 2],num[max(2*(i-n_num2),i), 1] - den[max(2*(i-n_den2),i), 1]];
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
          a := if abs(den[n_den2 + i, 1])>eps then -den[n_den2 + i, 1] else 0.0;

          if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
                // State space system in form (4)

            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]*(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/den[n_den2 + i, 1] else 1.0;
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1];
            d := k[i_k + i];
            dZero := false;
          elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
          // State space system in form (4)
            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]*(num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/den[n_den2 + i, 1] else 1.0;
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else num[max(1, n_num2 + i - i_d + 1), 1];
            d := k[i_k + i];
            dZero := false;
          else
           // State space system in form (5)
            c := if abs(den[n_den2 + i, 1])>eps then k[i_k + i]/den[n_den2 + i, 1] else k[i_k + i];
            b := if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1] else 1;
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
ss = ZerosAndPoles.Conversion.toStateSpace<b>toStateSpace</b>(zp)
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
        a2 + a1*p + p^2     a2     b2 + a1*b2/a2*p + b2/a2*p^2
G(p) = ----------------- = ---- * ----------------------------- = k * F(p),  k = a2/b2  (1)
        b2 + b1*p + p^2     b2           b2 + b1*p + p^2
</pre></blockquote>
<p>
for second order systems and
</p>
<blockquote><pre>
        a + p     a     b + b/a*p
G(p) = ------- = --- * ----------- = k * F(p),   k = a/b
        b + p     b      b + p
</pre></blockquote>
<p>
for first order systems respectively.
</p>
<p>
The complete system is now considered as the series connections of all
the single unit gain transfer functions and an overall gain k with
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
    "Package of functions to generate a ZerosAndPoles data record from imported data"
    extends Modelica.Icons.Package;
    import Modelica;

    encapsulated function fromFile
      "Generate a ZerosAndPoles data record by reading the polynomial coefficients or zeros and poles from a file"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DataDir;

      input String fileName=DataDir + "zp.mat"
        "Name of the zeros and poles data file"
        annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));

    protected
      Integer n1n2d1d2[4]=if ZerosAndPoles.Internal.checkRepresentation(
          fileName) then ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(
          fileName) else ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(
          fileName) annotation(__Dymola_allowForSize=true);
      Integer n1=n1n2d1d2[1] annotation(__Dymola_allowForSize=true);
      Integer n2=n1n2d1d2[2] annotation(__Dymola_allowForSize=true);
      Integer d1=n1n2d1d2[3] annotation(__Dymola_allowForSize=true);
      Integer d2=n1n2d1d2[4] annotation(__Dymola_allowForSize=true);
      Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2] annotation(__Dymola_allowForSize=true);
      Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4] annotation(__Dymola_allowForSize=true);
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
      //Whenever this function becomes operational the code must be
      // rewritten if fromFile_pc2 and fromFile_zp2 are in the 'constructor'

      zp := if ZerosAndPoles.Internal.checkRepresentation(fileName) then
        ZerosAndPoles.Internal.fromFile_zp(fileName) else
        ZerosAndPoles.Internal.fromFile_pc(fileName);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = ZerosAndPoles.Import.<b>fromFile</b>(fileName)
</pre></blockquote>

<h4>Description</h4>
<p>
Reads and loads a zeros-and-poles transfer function from a mat-file <tt>fileName</tt>. The file must contain either the set of variables n1, n2, d1, d2, and k with the associated first and second order polynomials or the variables p, z, and k with the poles and zeros, written in two column arrays with real and imaginary in the first and second column respectively. The variable k is the real gail in both cases.
</p>

<h4>Example</h4>
<blockquote><pre>
<b>algorithm</b>
  zp:=Modelica_LinearSystems2.ZerosAndPoles.Import.fromFile(\"zp.mat\", \"n\", \"d\");
//  zp = (p^2 + 2*p + 3)/(p + 2)/(p^2 + 2*p + 2)
</pre></blockquote>
</html>"));
    end fromFile;

    function fromModel
      "Generate a ZerosAndPoles data record from a state space representation resulted from linearization of a model"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

      input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel(translate=true)));
      input Real T_linearize = 0
        "Point in time of simulation to linearize the model";
      input String fileName = "dslin" "Name of the result file";

    protected
      String fileName2 = fileName + ".mat";
      Boolean OK1 = simulateModel(
            problem=modelName,
            startTime=0,
            stopTime=T_linearize);
      Boolean OK2 = importInitial("dsfinal.txt");
      Boolean OK3 = linearizeModel(
            problem=modelName,
            resultFile=fileName,
            startTime=T_linearize,
            stopTime=T_linearize + 1);
      Integer xuy[3] = ReadSystemDimension(fileName2, "ABCD");
      Integer nx = xuy[1];
      Integer nu = xuy[2];
      Integer ny = xuy[3];
      Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(
            fileName2,
            "ABCD",
            nx + ny,
            nx + nu);
      String xuyName[nx + nu + ny] = readStringMatrix(
            fileName2,
            "xuyName",
            nx + nu + ny);

      StateSpace result(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,nu],
        redeclare Real C[ny,nx],
        redeclare Real D[ny,nu]) "Model linearized at initial point";
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
<blockquote><pre>
zp = ZerosAndPoles.Import.<b>fromModel</b>(modelName, T_linearize, fileName)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate a matrix of ZerosAndPoles data records by linearization of a model
defined by modelName. The linearization is performed at time T_linearize of
the simulation. The system is genrated by using <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Import.fromFile\">StateSpace.Import.fromModel</a>
followed by a conversion from sate space to transfer function representation.
</p>

<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
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

  end Import;

  encapsulated package Internal
    "Package of internal material of record ZerosAndPoles (for advanced users only)"
    extends Modelica.Icons.InternalPackage;

    import Modelica;
    import Modelica_LinearSystems2;

  record ZerosAndPoles
      "Continuous zeros and poles description of a single input, single output system (data + operations)"
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

  end ZerosAndPoles;

    encapsulated function bandPassAlpha "Return alpha for band pass"
      import Modelica;

      input Real a "Coefficient of p^1";
      input Real b "Coefficient of p^0";
      input Modelica.SIunits.AngularVelocity w "Bandwidth angular frequency";
      output Real alpha "Alpha factor to build up band pass";

    protected
      Real alpha_min;
      Real alpha_max;
      Real z_min;
      Real z_max;
      Real z;

      function residue "Residue of non-linear equation"
        input Real a;
        input Real b;
        input Real w;
        input Real z;
        output Real res;
      algorithm
        res := z^2 + (a*w*z/(1+z))^2 - (2+b*w^2)*z + 1;
      end residue;

    function solveOneNonlinearEquation
        "Solve f(u) = 0; f(u_min) and f(u_max) must have different signs"
        import Modelica.Utilities.Streams.error;

      input Real aa;
      input Real bb;
      input Real ww;
      input Real u_min "Lower bound of search intervall";
      input Real u_max "Upper bound of search intervall";
      input Real tolerance=100*Modelica.Constants.eps
          "Relative tolerance of solution u";
      output Real u "Value of independent variable so that f(u) = 0";

      protected
      constant Real eps=Modelica.Constants.eps "machine epsilon";
      Real a=u_min "Current best minimum interval value";
      Real b=u_max "Current best maximum interval value";
      Real c "Intermediate point a <= c <= b";
      Real d;
      Real e "b - a";
      Real m;
      Real s;
      Real p;
      Real q;
      Real r;
      Real tol;
      Real fa "= f(a)";
      Real fb "= f(b)";
      Real fc;
      Boolean found=false;
    algorithm
      // Check that f(u_min) and f(u_max) have different sign
      fa := residue(aa,bb,ww,u_min);
      fb := residue(aa,bb,ww,u_max);
      fc := fb;
      if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
        error(
          "The arguments u_min and u_max to solveOneNonlinearEquation(..)\n" +
          "do not bracket the root of the single non-linear equation:\n" +
          "  u_min  = " + String(u_min) + "\n" + "  u_max  = " + String(u_max)
           + "\n" + "  fa = f(u_min) = " + String(fa) + "\n" +
          "  fb = f(u_max) = " + String(fb) + "\n" +
          "fa and fb must have opposite sign which is not the case");
      end if;

      // Initialize variables
      c := a;
      fc := fa;
      e := b - a;
      d := e;

      // Search loop
      while not found loop
        if abs(fc) < abs(fb) then
          a := b;
          b := c;
          c := a;
          fa := fb;
          fb := fc;
          fc := fa;
        end if;

        tol := 2*eps*abs(b) + tolerance;
        m := (c - b)/2;

        if abs(m) <= tol or fb == 0.0 then
          // root found (interval is small enough)
          found := true;
          u := b;
        else
          // Determine if a bisection is needed
          if abs(e) < tol or abs(fa) <= abs(fb) then
            e := m;
            d := e;
          else
            s := fb/fa;
            if a == c then
              // linear interpolation
              p := 2*m*s;
              q := 1 - s;
            else
              // inverse quadratic interpolation
              q := fa/fc;
              r := fb/fc;
              p := s*(2*m*q*(q - r) - (b - a)*(r - 1));
              q := (q - 1)*(r - 1)*(s - 1);
            end if;

            if p > 0 then
              q := -q;
            else
              p := -p;
            end if;

            s := e;
            e := d;
            if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q) then
              // interpolation successful
              d := p/q;
            else
              // use bi-section
              e := m;
              d := e;
            end if;
          end if;

          // Best guess value is defined as "a"
          a := b;
          fa := fb;
          b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
          fb := residue(aa,bb,ww,b);

          if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
            // initialize variables
            c := a;
            fc := fa;
            e := b - a;
            d := e;
          end if;
        end if;
      end while;

      annotation (Documentation(info="<html>

<p>
This function determines the solution of <b>one non-linear algebraic equation</b> &quot;y=f(u)&quot;
in <b>one unknown</b> &quot;u&quot; in a reliable way. It is one of the best numerical
algorithms for this purpose. As input, the nonlinear function f(u)
has to be given, as well as an interval u_min, u_max that
contains the solution, i.e., &quot;f(u_min)&quot; and &quot;f(u_max)&quot; must
have a different sign. If possible, a smaller interval is computed by
inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last 3 points and computing the zero). If this fails,
bisection is used, which always reduces the interval by a factor of 2.
The inverse quadratic interpolation method has superlinear convergence.
This is roughly the same convergence rate as a globally convergent Newton
method, but without the need to compute derivatives of the non-linear
function. The solver function is a direct mapping of the Algol 60 procedure
&quot;zero&quot; to Modelica, from:
</p>

<dl>
<dt> Brent R.P.:</dt>
<dd> <b>Algorithms for Minimization without derivatives</b>.
     Prentice Hall, 1973, pp. 58-59.</dd>
</dl>
</html>"));
    end solveOneNonlinearEquation;

    algorithm
      assert( a^2/4 - b <= 0,  "Band pass transformation cannot be computed");
      z :=solveOneNonlinearEquation(a, b, w, 0, 1);
      alpha := sqrt(z);

      annotation (Documentation(info="<html>
<p>
A band pass with bandwidth &quot;w&quot; is determined from a low pass
</p>

<blockquote><pre>
  1/(p^2 + a*p + b)
</pre></blockquote>

<p>
with the transformation
</p>

<blockquote><pre>
  new(p) = (p + 1/p)/w
</pre></blockquote>

<p>
This results in the following derivation:
</p>

<blockquote><pre>
  1/(p^2 + a*p + b) -> 1/( (p+1/p)^2/w^2 + a*(p + 1/p)/w + b )
                     = 1 / ( p^2 + 1/p^2 + 2)/w^2 + (p + 1/p)*a/w + b )
                     = w^2*p^2 / (p^4 + 2*p^2 + 1 + (p^3 + p)a*w + b*w^2*p^2)
                     = w^2*p^2 / (p^4 + a*w*p^3 + (2+b*w^2)*p^2 + a*w*p + 1)
</pre></blockquote>

<p>
This 4th order transfer function shall be split in to two transfer functions of order 2 each
for numerical reasons. With the following formulation, the fourth order
polynomial can be represented (with the unknowns &quot;c&quot; and &quot;alpha&quot;):
</p>

<blockquote><pre>
  g(p) = w^2*p^2 / ( (p*alpha)^2 + c*(p*alpha) + 1) * (p/alpha)^2 + c*(p/alpha) + 1)
       = w^2*p^2 / ( p^4 + c*(alpha + 1/alpha)*p^3 + (alpha^2 + 1/alpha^2 + c^2)*p^2
                                                   + c*(alpha + 1/alpha)*p + 1 )
</pre></blockquote>

<p>
Comparison of coefficients:
</p>

<blockquote><pre>
  c*(alpha + 1/alpha) = a*w           -> c = a*w / (alpha + 1/alpha)
  alpha^2 + 1/alpha^2 + c^2 = 2+b*w^2 -> equation to determine alpha

  alpha^4 + 1 + a^2*w^2*alpha^4/(1+alpha^2)^2 = (2+b*w^2)*alpha^2
    or z = alpha^2
  z^2 + a^2*w^2*z^2/(1+z)^2 - (2+b*w^2)*z + 1 = 0
</pre></blockquote>

<p>
Therefore the last equation has to be solved for &quot;z&quot; (basically, this means to compute
a real zero of a fourth order polynomal):
</p>

<blockquote><pre>
   solve: 0 = f(z)  = z^2 + a^2*w^2*z^2/(1+z)^2 - (2+b*w^2)*z + 1  for &quot;z&quot;
              f(0)  = 1  &gt; 0
              f(1)  = 1 + a^2*w^2/4 - (2+b*w^2) + 1
                    = (a^2/4 - b)*w^2  // must be &lt; 0
</pre></blockquote>

<p>
This function computes the solution of this equation and returns &quot;alpha = z^2&quot; and &quot;c&quot;;
</p>
</html>"));
    end bandPassAlpha;

  encapsulated function baseFilter
      "Generate a ZerosAndPoles transfer function from a base filter description (= low pass filter with w_cut = 1 rad/s)"

      import Modelica;
      import Modelica.Math;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Utilities.Types;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input Modelica_LinearSystems2.Utilities.Types.AnalogFilter analogFilter=Modelica_LinearSystems2.Utilities.Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
    input Integer order(min=1) = 2 "Order of filter";
    input Real A_ripple(unit="dB") = 0.5
        "Pass band ripple (only for Chebyshev filter)";
    input Boolean normalized=true
        "True, if amplitude at f_cut = -3db, otherwise unmodified filter";

    output ZerosAndPoles filter(
      redeclare Real n1[0],
      redeclare Real n2[0,2],
      redeclare Real d1[if analogFilter == Types.AnalogFilter.CriticalDamping then
              order else mod(order, 2)],
      redeclare Real d2[if analogFilter == Types.AnalogFilter.CriticalDamping then
              0 else integer(order/2),2]) "Filter transfer function";
    protected
    Integer n_den1=size(filter.d1, 1);
    Integer n_den2=size(filter.d2, 1);
    Integer n_den=n_den1 + 2*n_den2;
    Real pi=Modelica.Constants.pi;
    Boolean evenOrder=mod(order, 2) == 0
        "True, if even filter order (otherwise uneven)";
    Real alpha=1.0 "Frequency correction factor";
    Real alpha2 "= alpha*alpha";
    Real epsilon "Ripple size";
    Real fac "arsinh(epsilon)";
    Real den1[n_den1]
        "[p] coefficients of denominator first order polynomials (a*p + 1)";
    Real den2[n_den2,2]
        "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1)";
    Real aux;
    Real k;
  algorithm
    /* Compute filter coefficients of prototype low pass filter. */
    if analogFilter == Types.AnalogFilter.CriticalDamping then
      if normalized then
        alpha := sqrt(2^(1/order) - 1);
        // alpha := sqrt(10^(3/10/order)-1);
      else
        alpha := 1;
      end if;
      for i in 1:n_den1 loop
        den1[i] := alpha;
      end for;

    elseif analogFilter == Types.AnalogFilter.Bessel then
      (den1,den2,alpha) := ZerosAndPoles.Internal.BesselCoefficients(order);
      if not normalized then
        alpha2 := alpha*alpha;
        for i in 1:n_den2 loop
          den2[i, 1] := den2[i, 1]*alpha2;
          den2[i, 2] := den2[i, 2]*alpha;
        end for;
        if not evenOrder then
          den1[1] := den1[1]*alpha;
        end if;
      end if;

    elseif analogFilter == Types.AnalogFilter.Butterworth then
       // Original filter is already normalized
      for i in 1:n_den2 loop
        den2[i, 1] := 1.0;
        den2[i, 2] := -2*cos(pi*(0.5 + (i - 0.5)/order));
      end for;
      if not evenOrder then
        den1[1] := 1.0;
      end if;

      /* Transformation of filter transfer function with "new(p) = alpha*p"
       in order that the filter transfer function has an amplitude of
       -3 db at the cutoff frequency
    */
      /*
    if normalized then
      alpha := ZerosAndPoles.Internal.normalizationFactor(den1, den2);
      alpha2 := alpha*alpha;
      for i in 1:n_den2 loop
        den2[i, 1] := den2[i, 1]*alpha2;
        den2[i, 2] := den2[i, 2]*alpha;
      end for;
      if not evenOrder then
        den1[1] := den1[1]*alpha;
      end if;
    end if;
    */

    elseif analogFilter == Types.AnalogFilter.Chebyshev then
      epsilon := sqrt(10^(A_ripple/10) - 1);
      fac := Math.asinh(1/epsilon)/order;

      if evenOrder then
         for i in 1:n_den2 loop
            den2[i,1] :=1/(cosh(fac)^2 - cos((2*i - 1)*pi/(2*order))^2);
            den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos((2*i - 1)*pi/(2*order));
         end for;
      else
         den1[1] := 1/sinh(fac);
         for i in 1:n_den2 loop
            den2[i,1] :=1/(cosh(fac)^2 - cos(i*pi/order)^2);
            den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos(i*pi/order);
         end for;
      end if;

       /* Transformation of filter transfer function with "new(p) = alpha*p"
        in order that the filter transfer function has an amplitude of
        -3 db at the cutoff frequency
     */
      if normalized then
        alpha := ZerosAndPoles.Internal.normalizationFactor(den1, den2);
        alpha2 := alpha*alpha;
        for i in 1:n_den2 loop
          den2[i, 1] := den2[i, 1]*alpha2;
          den2[i, 2] := den2[i, 2]*alpha;
        end for;
        if not evenOrder then
          den1[1] := den1[1]*alpha;
        end if;
      end if;

    else
      Streams.error("analogFilter (= " + String(analogFilter) +
        ") is not supported");
    end if;

    // Determine normalized denominator polynomials with highest power of p equal to one
    (filter.d1,filter.d2,k) := ZerosAndPoles.Internal.filterToNormalized(den1, den2);
    filter.k := 1.0/k;

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = <b>baseFilter</b>(analogFilter, order, A_ripple, normalized);
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles transfer function
description of low pass filters with a cut-off angular frequency
of one rad/s and a gain of one. Filters returned by this function
are the starting point to construct other filters by transformation
of the filter transfer function:
</p>

<blockquote><pre>
zp(p) = 1 / ( product( p + a[i] ) * product(p^2 + b[i]*p + a[i]) )
</pre></blockquote>

<p>
using the following rules:
</p>

<table border=1 cellspacing=0 cellpadding=2>
<tr><td><i>Desired filter</i></td>
    <td><i>Transformation</i></td>
    </tr>

<tr><td> High pass filter </td>
    <td> replace &quot;p&quot; by &quot;1/p&quot; </td>
    </tr>

<tr><td> Band pass filter </td>
    <td> replace &quot;p&quot; by &quot;(p + 1/p)/w_band&quot;<br>
         (w_band = (f_max - f_min)/sqrt(f_min*f_max))</td>
    </tr>

<tr><td> Stop pass filter  </td>
    <td> replace &quot;p&quot; by &quot;w_band/(p + 1/p)&quot;<br>
         (w_band = (f_max - f_min)/sqrt(f_min*f_max))</td>
    </tr>

<tr><td> Filter with cut-off angular frequency w_cut </td>
    <td> replace &quot;p&quot; by &quot;p/w_cut&quot; </td>
    </tr>
</table>
<p>
For more details see also
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Literature\">[Tietze2002]</a>, pp. 815-852.
</p>

<h4>Example</h4>
<blockquote><pre>
  // Generate a Butterworth high pass base filter of order 3
  <b>import</b> ZP = Modelica_LinearSystems2.ZerosAndPoles;
  <b>import</b> Modelica_LinearSystems2.Types;

  ZP zp_filter;
<b>algorithm</b>
  zp_filter = ZP.Internal.baseFilter(Types.AnalogFilter.Butterworth, order = 3);

 // zp_filter = 1 /  ( (p + 1)*(p^2 + p + 1) )
</pre></blockquote>
</html>"));
  end baseFilter;

    function BesselCoefficients
      "Return coefficients of normalized low pass Bessel filter (= gain at cut-off frequency 1 rad/s is decreased 3dB"

      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      input Integer order "Order of filter in the range 1..41";
      output Real c1[mod(order, 2)]
        "[p] coefficients of Bessel denominator polynomials (a*p + 1)";
      output Real c2[integer(order/2),2]
        "[p^2, p] coefficients of Bessel denominator polynomials (b2*p^2 + b1*p + 1)";
      output Real alpha "Normalization factor";
    algorithm
      if order == 1 then
        alpha := 1.002377293007601;
        c1[1] := 0.9976283451109835;
      elseif order == 2 then
        alpha := 0.7356641785819585;
        c2[1, 1] := 0.6159132201783791;
        c2[1, 2] := 1.359315879600889;
      elseif order == 3 then
        alpha := 0.5704770156982642;
        c1[1] := 0.7548574865985343;
        c2[1, 1] := 0.4756958028827457;
        c2[1, 2] := 0.9980615136104388;
      elseif order == 4 then
        alpha := 0.4737978580281427;
        c2[1, 1] := 0.4873729247240677;
        c2[1, 2] := 1.337564170455762;
        c2[2, 1] := 0.3877724315741958;
        c2[2, 2] := 0.7730405590839861;
      elseif order == 5 then
        alpha := 0.4126226974763408;
        c1[1] := 0.6645723262620757;
        c2[1, 1] := 0.4115231900614016;
        c2[1, 2] := 1.138349926728708;
        c2[2, 1] := 0.3234938702877912;
        c2[2, 2] := 0.6205992985771313;
      elseif order == 6 then
        alpha := 0.3705098000736233;
        c2[1, 1] := 0.3874508649098960;
        c2[1, 2] := 1.219740879520741;
        c2[2, 1] := 0.3493298843155746;
        c2[2, 2] := 0.9670265529381365;
        c2[3, 1] := 0.2747419229514599;
        c2[3, 2] := 0.5122165075105700;
      elseif order == 7 then
        alpha := 0.3393452623586350;
        c1[1] := 0.5927147125821412;
        c2[1, 1] := 0.3383379423919174;
        c2[1, 2] := 1.092630816438030;
        c2[2, 1] := 0.3001025788696046;
        c2[2, 2] := 0.8289928256598656;
        c2[3, 1] := 0.2372867471539579;
        c2[3, 2] := 0.4325128641920154;
      elseif order == 8 then
        alpha := 0.3150267393795002;
        c2[1, 1] := 0.3151115975207653;
        c2[1, 2] := 1.109403015460190;
        c2[2, 1] := 0.2969344839572762;
        c2[2, 2] := 0.9737455812222699;
        c2[3, 1] := 0.2612545921889538;
        c2[3, 2] := 0.7190394712068573;
        c2[4, 1] := 0.2080523342974281;
        c2[4, 2] := 0.3721456473047434;
      elseif order == 9 then
        alpha := 0.2953310177184124;
        c1[1] := 0.5377196679501422;
        c2[1, 1] := 0.2824689124281034;
        c2[1, 2] := 1.022646191567475;
        c2[2, 1] := 0.2626824161383468;
        c2[2, 2] := 0.8695626454762596;
        c2[3, 1] := 0.2302781917677917;
        c2[3, 2] := 0.6309047553448520;
        c2[4, 1] := 0.1847991729757028;
        c2[4, 2] := 0.3251978031287202;
      elseif order == 10 then
        alpha := 0.2789426890619463;
        c2[1, 1] := 0.2640769908255582;
        c2[1, 2] := 1.019788132875305;
        c2[2, 1] := 0.2540802639216947;
        c2[2, 2] := 0.9377020417760623;
        c2[3, 1] := 0.2343577229427963;
        c2[3, 2] := 0.7802229808216112;
        c2[4, 1] := 0.2052193139338624;
        c2[4, 2] := 0.5594176813008133;
        c2[5, 1] := 0.1659546953748916;
        c2[5, 2] := 0.2878349616233292;
      elseif order == 11 then
        alpha := 0.2650227766037203;
        c1[1] := 0.4950265498954191;
        c2[1, 1] := 0.2411858478546218;
        c2[1, 2] := 0.9567800996387417;
        c2[2, 1] := 0.2296849355380925;
        c2[2, 2] := 0.8592523717113126;
        c2[3, 1] := 0.2107851705677406;
        c2[3, 2] := 0.7040216048898129;
        c2[4, 1] := 0.1846461385164021;
        c2[4, 2] := 0.5006729207276717;
        c2[5, 1] := 0.1504217970817433;
        c2[5, 2] := 0.2575070491320295;
      elseif order == 12 then
        alpha := 0.2530051198547209;
        c2[1, 1] := 0.2268294941204543;
        c2[1, 2] := 0.9473116570034053;
        c2[2, 1] := 0.2207657387793729;
        c2[2, 2] := 0.8933728946287606;
        c2[3, 1] := 0.2087600700376653;
        c2[3, 2] := 0.7886236252756229;
        c2[4, 1] := 0.1909959101492760;
        c2[4, 2] := 0.6389263649257017;
        c2[5, 1] := 0.1675208146048472;
        c2[5, 2] := 0.4517847275162215;
        c2[6, 1] := 0.1374257286372761;
        c2[6, 2] := 0.2324699157474680;
      elseif order == 13 then
        alpha := 0.2424910397561007;
        c1[1] := 0.4608848369928040;
        c2[1, 1] := 0.2099813050274780;
        c2[1, 2] := 0.8992478823790660;
        c2[2, 1] := 0.2027250423101359;
        c2[2, 2] := 0.8328117484224146;
        c2[3, 1] := 0.1907635894058731;
        c2[3, 2] := 0.7257379204691213;
        c2[4, 1] := 0.1742280397887686;
        c2[4, 2] := 0.5830640944868014;
        c2[5, 1] := 0.1530858190490478;
        c2[5, 2] := 0.4106192089751885;
        c2[6, 1] := 0.1264090712880446;
        c2[6, 2] := 0.2114980230156001;
      elseif order == 14 then
        alpha := 0.2331902368695848;
        c2[1, 1] := 0.1986162311411235;
        c2[1, 2] := 0.8876961808055535;
        c2[2, 1] := 0.1946683341271615;
        c2[2, 2] := 0.8500754229171967;
        c2[3, 1] := 0.1868331332895056;
        c2[3, 2] := 0.7764629313723603;
        c2[4, 1] := 0.1752118757862992;
        c2[4, 2] := 0.6699720402924552;
        c2[5, 1] := 0.1598906457908402;
        c2[5, 2] := 0.5348446712848934;
        c2[6, 1] := 0.1407810153019944;
        c2[6, 2] := 0.3755841316563539;
        c2[7, 1] := 0.1169627966707339;
        c2[7, 2] := 0.1937088226304455;
      elseif order == 15 then
        alpha := 0.2248854870552422;
        c1[1] := 0.4328492272335646;
        c2[1, 1] := 0.1857292591004588;
        c2[1, 2] := 0.8496337061962563;
        c2[2, 1] := 0.1808644178280136;
        c2[2, 2] := 0.8020517898136011;
        c2[3, 1] := 0.1728264404199081;
        c2[3, 2] := 0.7247449729331105;
        c2[4, 1] := 0.1616970125901954;
        c2[4, 2] := 0.6205369315943097;
        c2[5, 1] := 0.1475257264578426;
        c2[5, 2] := 0.4929612162355906;
        c2[6, 1] := 0.1301861023357119;
        c2[6, 2] := 0.3454770708040735;
        c2[7, 1] := 0.1087810777120188;
        c2[7, 2] := 0.1784526655428406;
      elseif order == 16 then
        alpha := 0.2174105053474761;
        c2[1, 1] := 0.1765637967473151;
        c2[1, 2] := 0.8377453068635511;
        c2[2, 1] := 0.1738525357503125;
        c2[2, 2] := 0.8102988957433199;
        c2[3, 1] := 0.1684627004613343;
        c2[3, 2] := 0.7563265923413258;
        c2[4, 1] := 0.1604519074815815;
        c2[4, 2] := 0.6776082294687619;
        c2[5, 1] := 0.1498828607802206;
        c2[5, 2] := 0.5766417034027680;
        c2[6, 1] := 0.1367764717792823;
        c2[6, 2] := 0.4563528264410489;
        c2[7, 1] := 0.1209810465419295;
        c2[7, 2] := 0.3193782657322374;
        c2[8, 1] := 0.1016312648007554;
        c2[8, 2] := 0.1652419227369036;
      elseif order == 17 then
        alpha := 0.2106355148193306;
        c1[1] := 0.4093223608497299;
        c2[1, 1] := 0.1664014345826274;
        c2[1, 2] := 0.8067173752345952;
        c2[2, 1] := 0.1629839591538256;
        c2[2, 2] := 0.7712924931447541;
        c2[3, 1] := 0.1573277802512491;
        c2[3, 2] := 0.7134213666303411;
        c2[4, 1] := 0.1494828185148637;
        c2[4, 2] := 0.6347841731714884;
        c2[5, 1] := 0.1394948812681826;
        c2[5, 2] := 0.5375594414619047;
        c2[6, 1] := 0.1273627583380806;
        c2[6, 2] := 0.4241608926375478;
        c2[7, 1] := 0.1129187258461290;
        c2[7, 2] := 0.2965752009703245;
        c2[8, 1] := 0.9533357359908857e-1;
        c2[8, 2] := 0.1537041700889585;
      elseif order == 18 then
        alpha := 0.2044575288651841;
        c2[1, 1] := 0.1588768571976356;
        c2[1, 2] := 0.7951914263212913;
        c2[2, 1] := 0.1569357024981854;
        c2[2, 2] := 0.7744529690772538;
        c2[3, 1] := 0.1530722206358810;
        c2[3, 2] := 0.7335304425992080;
        c2[4, 1] := 0.1473206710524167;
        c2[4, 2] := 0.6735038935387268;
        c2[5, 1] := 0.1397225420331520;
        c2[5, 2] := 0.5959151542621590;
        c2[6, 1] := 0.1303092459809849;
        c2[6, 2] := 0.5026483447894845;
        c2[7, 1] := 0.1190627367060072;
        c2[7, 2] := 0.3956893824587150;
        c2[8, 1] := 0.1058058030798994;
        c2[8, 2] := 0.2765091830730650;
        c2[9, 1] := 0.8974708108800873e-1;
        c2[9, 2] := 0.1435505288284833;
      elseif order == 19 then
        alpha := 0.1987936248083529;
        c1[1] := 0.3892259966869526;
        c2[1, 1] := 0.1506640012172225;
        c2[1, 2] := 0.7693121733774260;
        c2[2, 1] := 0.1481728062796673;
        c2[2, 2] := 0.7421133586741549;
        c2[3, 1] := 0.1440444668388838;
        c2[3, 2] := 0.6975075386214800;
        c2[4, 1] := 0.1383101628540374;
        c2[4, 2] := 0.6365464378910025;
        c2[5, 1] := 0.1310032283190998;
        c2[5, 2] := 0.5606211948462122;
        c2[6, 1] := 0.1221431166405330;
        c2[6, 2] := 0.4713530424221445;
        c2[7, 1] := 0.1116991161103884;
        c2[7, 2] := 0.3703717538617073;
        c2[8, 1] := 0.9948917351196349e-1;
        c2[8, 2] := 0.2587371155559744;
        c2[9, 1] := 0.8475989238107367e-1;
        c2[9, 2] := 0.1345537894555993;
      elseif order == 20 then
        alpha := 0.1935761760416219;
        c2[1, 1] := 0.1443871348337404;
        c2[1, 2] := 0.7584165598446141;
        c2[2, 1] := 0.1429501891353184;
        c2[2, 2] := 0.7423000962318863;
        c2[3, 1] := 0.1400877384920004;
        c2[3, 2] := 0.7104185332215555;
        c2[4, 1] := 0.1358210369491446;
        c2[4, 2] := 0.6634599783272630;
        c2[5, 1] := 0.1301773703034290;
        c2[5, 2] := 0.6024175491895959;
        c2[6, 1] := 0.1231826501439148;
        c2[6, 2] := 0.5285332736326852;
        c2[7, 1] := 0.1148465498575254;
        c2[7, 2] := 0.4431977385498628;
        c2[8, 1] := 0.1051289462376788;
        c2[8, 2] := 0.3477444062821162;
        c2[9, 1] := 0.9384622797485121e-1;
        c2[9, 2] := 0.2429038300327729;
        c2[10, 1] := 0.8028211612831444e-1;
        c2[10, 2] := 0.1265329974009533;
      elseif order == 21 then
        alpha := 0.1887494014766075;
        c1[1] := 0.3718070668941645;
        c2[1, 1] := 0.1376151928386445;
        c2[1, 2] := 0.7364290859445481;
        c2[2, 1] := 0.1357438914390695;
        c2[2, 2] := 0.7150167318935022;
        c2[3, 1] := 0.1326398453462415;
        c2[3, 2] := 0.6798001808470175;
        c2[4, 1] := 0.1283231214897678;
        c2[4, 2] := 0.6314663440439816;
        c2[5, 1] := 0.1228169159777534;
        c2[5, 2] := 0.5709353626166905;
        c2[6, 1] := 0.1161406100773184;
        c2[6, 2] := 0.4993087153571335;
        c2[7, 1] := 0.1082959649233524;
        c2[7, 2] := 0.4177766148584385;
        c2[8, 1] := 0.9923596957485723e-1;
        c2[8, 2] := 0.3274257287232124;
        c2[9, 1] := 0.8877776108724853e-1;
        c2[9, 2] := 0.2287218166767916;
        c2[10, 1] := 0.7624076527736326e-1;
        c2[10, 2] := 0.1193423971506988;
      elseif order == 22 then
        alpha := 0.1842668221199706;
        c2[1, 1] := 0.1323053462701543;
        c2[1, 2] := 0.7262446126765204;
        c2[2, 1] := 0.1312121721769772;
        c2[2, 2] := 0.7134286088450949;
        c2[3, 1] := 0.1290330911166814;
        c2[3, 2] := 0.6880287870435514;
        c2[4, 1] := 0.1257817990372067;
        c2[4, 2] := 0.6505015800059301;
        c2[5, 1] := 0.1214765261983008;
        c2[5, 2] := 0.6015107185211451;
        c2[6, 1] := 0.1161365140967959;
        c2[6, 2] := 0.5418983553698413;
        c2[7, 1] := 0.1097755171533100;
        c2[7, 2] := 0.4726370779831614;
        c2[8, 1] := 0.1023889478519956;
        c2[8, 2] := 0.3947439506537486;
        c2[9, 1] := 0.9392485861253800e-1;
        c2[9, 2] := 0.3090996703083202;
        c2[10, 1] := 0.8420273775456455e-1;
        c2[10, 2] := 0.2159561978556017;
        c2[11, 1] := 0.7257600023938262e-1;
        c2[11, 2] := 0.1128633732721116;
      elseif order == 23 then
        alpha := 0.1800893554453722;
        c1[1] := 0.3565232673929280;
        c2[1, 1] := 0.1266275171652706;
        c2[1, 2] := 0.7072778066734162;
        c2[2, 1] := 0.1251865227648538;
        c2[2, 2] := 0.6900676345785905;
        c2[3, 1] := 0.1227944815236645;
        c2[3, 2] := 0.6617011100576023;
        c2[4, 1] := 0.1194647013077667;
        c2[4, 2] := 0.6226432315773119;
        c2[5, 1] := 0.1152132989252356;
        c2[5, 2] := 0.5735222810625359;
        c2[6, 1] := 0.1100558598478487;
        c2[6, 2] := 0.5151027978024605;
        c2[7, 1] := 0.1040013558214886;
        c2[7, 2] := 0.4482410942032739;
        c2[8, 1] := 0.9704014176512626e-1;
        c2[8, 2] := 0.3738049984631116;
        c2[9, 1] := 0.8911683905758054e-1;
        c2[9, 2] := 0.2925028692588410;
        c2[10, 1] := 0.8005438265072295e-1;
        c2[10, 2] := 0.2044134600278901;
        c2[11, 1] := 0.6923832296800832e-1;
        c2[11, 2] := 0.1069984887283394;
      elseif order == 24 then
        alpha := 0.1761838665838427;
        c2[1, 1] := 0.1220804912720132;
        c2[1, 2] := 0.6978026874156063;
        c2[2, 1] := 0.1212296762358897;
        c2[2, 2] := 0.6874139794926736;
        c2[3, 1] := 0.1195328372961027;
        c2[3, 2] := 0.6667954259551859;
        c2[4, 1] := 0.1169990987333593;
        c2[4, 2] := 0.6362602049901176;
        c2[5, 1] := 0.1136409040480130;
        c2[5, 2] := 0.5962662188435553;
        c2[6, 1] := 0.1094722001757955;
        c2[6, 2] := 0.5474001634109253;
        c2[7, 1] := 0.1045052832229087;
        c2[7, 2] := 0.4903523180249535;
        c2[8, 1] := 0.9874509806025907e-1;
        c2[8, 2] := 0.4258751523524645;
        c2[9, 1] := 0.9217799943472177e-1;
        c2[9, 2] := 0.3547079765396403;
        c2[10, 1] := 0.8474633796250476e-1;
        c2[10, 2] := 0.2774145482392767;
        c2[11, 1] := 0.7627722381240495e-1;
        c2[11, 2] := 0.1939329108084139;
        c2[12, 1] := 0.6618645465422745e-1;
        c2[12, 2] := 0.1016670147947242;
      elseif order == 25 then
        alpha := 0.1725220521949266;
        c1[1] := 0.3429735385896000;
        c2[1, 1] := 0.1172525033170618;
        c2[1, 2] := 0.6812327932576614;
        c2[2, 1] := 0.1161194585333535;
        c2[2, 2] := 0.6671566071153211;
        c2[3, 1] := 0.1142375145794466;
        c2[3, 2] := 0.6439167855053158;
        c2[4, 1] := 0.1116157454252308;
        c2[4, 2] := 0.6118378416180135;
        c2[5, 1] := 0.1082654809459177;
        c2[5, 2] := 0.5713609763370088;
        c2[6, 1] := 0.1041985674230918;
        c2[6, 2] := 0.5230289949762722;
        c2[7, 1] := 0.9942439308123559e-1;
        c2[7, 2] := 0.4674627926041906;
        c2[8, 1] := 0.9394453593830893e-1;
        c2[8, 2] := 0.4053226688298811;
        c2[9, 1] := 0.8774221237222533e-1;
        c2[9, 2] := 0.3372372276379071;
        c2[10, 1] := 0.8075839512216483e-1;
        c2[10, 2] := 0.2636485508005428;
        c2[11, 1] := 0.7282483286646764e-1;
        c2[11, 2] := 0.1843801345273085;
        c2[12, 1] := 0.6338571166846652e-1;
        c2[12, 2] := 0.9680153764737715e-1;
      elseif order == 26 then
        alpha := 0.1690795702796737;
        c2[1, 1] := 0.1133168695796030;
        c2[1, 2] := 0.6724297955493932;
        c2[2, 1] := 0.1126417845769961;
        c2[2, 2] := 0.6638709519790540;
        c2[3, 1] := 0.1112948749545606;
        c2[3, 2] := 0.6468652038763624;
        c2[4, 1] := 0.1092823986944244;
        c2[4, 2] := 0.6216337070799265;
        c2[5, 1] := 0.1066130386697976;
        c2[5, 2] := 0.5885011413992190;
        c2[6, 1] := 0.1032969057045413;
        c2[6, 2] := 0.5478864278297548;
        c2[7, 1] := 0.9934388184210715e-1;
        c2[7, 2] := 0.5002885306054287;
        c2[8, 1] := 0.9476081523436283e-1;
        c2[8, 2] := 0.4462644847551711;
        c2[9, 1] := 0.8954648464575577e-1;
        c2[9, 2] := 0.3863930785049522;
        c2[10, 1] := 0.8368166847159917e-1;
        c2[10, 2] := 0.3212074592527143;
        c2[11, 1] := 0.7710664731701103e-1;
        c2[11, 2] := 0.2510470347119383;
        c2[12, 1] := 0.6965807988411425e-1;
        c2[12, 2] := 0.1756419294111342;
        c2[13, 1] := 0.6080674930548766e-1;
        c2[13, 2] := 0.9234535279274277e-1;
      elseif order == 27 then
        alpha := 0.1658353543067995;
        c1[1] := 0.3308543720638957;
        c2[1, 1] := 0.1091618578712746;
        c2[1, 2] := 0.6577977071169651;
        c2[2, 1] := 0.1082549561495043;
        c2[2, 2] := 0.6461121666520275;
        c2[3, 1] := 0.1067479247890451;
        c2[3, 2] := 0.6267937760991321;
        c2[4, 1] := 0.1046471079537577;
        c2[4, 2] := 0.6000750116745808;
        c2[5, 1] := 0.1019605976654259;
        c2[5, 2] := 0.5662734183049320;
        c2[6, 1] := 0.9869726954433709e-1;
        c2[6, 2] := 0.5257827234948534;
        c2[7, 1] := 0.9486520934132483e-1;
        c2[7, 2] := 0.4790595019077763;
        c2[8, 1] := 0.9046906518775348e-1;
        c2[8, 2] := 0.4266025862147336;
        c2[9, 1] := 0.8550529998276152e-1;
        c2[9, 2] := 0.3689188223512328;
        c2[10, 1] := 0.7995282239306020e-1;
        c2[10, 2] := 0.3064589322702932;
        c2[11, 1] := 0.7375174596252882e-1;
        c2[11, 2] := 0.2394754504667310;
        c2[12, 1] := 0.6674377263329041e-1;
        c2[12, 2] := 0.1676223546666024;
        c2[13, 1] := 0.5842458027529246e-1;
        c2[13, 2] := 0.8825044329219431e-1;
      elseif order == 28 then
        alpha := 0.1627710671942929;
        c2[1, 1] := 0.1057232656113488;
        c2[1, 2] := 0.6496161226860832;
        c2[2, 1] := 0.1051786825724864;
        c2[2, 2] := 0.6424661279909941;
        c2[3, 1] := 0.1040917964935006;
        c2[3, 2] := 0.6282470268918791;
        c2[4, 1] := 0.1024670101953951;
        c2[4, 2] := 0.6071189030701136;
        c2[5, 1] := 0.1003105109519892;
        c2[5, 2] := 0.5793175191747016;
        c2[6, 1] := 0.9762969425430802e-1;
        c2[6, 2] := 0.5451486608855443;
        c2[7, 1] := 0.9443223803058400e-1;
        c2[7, 2] := 0.5049796971628137;
        c2[8, 1] := 0.9072460982036488e-1;
        c2[8, 2] := 0.4592270546572523;
        c2[9, 1] := 0.8650956423253280e-1;
        c2[9, 2] := 0.4083368605952977;
        c2[10, 1] := 0.8178165740374893e-1;
        c2[10, 2] := 0.3527525188880655;
        c2[11, 1] := 0.7651838885868020e-1;
        c2[11, 2] := 0.2928534570013572;
        c2[12, 1] := 0.7066010532447490e-1;
        c2[12, 2] := 0.2288185204390681;
        c2[13, 1] := 0.6405358596145789e-1;
        c2[13, 2] := 0.1602396172588190;
        c2[14, 1] := 0.5621780070227172e-1;
        c2[14, 2] := 0.8447589564915071e-1;
      elseif order == 29 then
        alpha := 0.1598706626277596;
        c1[1] := 0.3199314513011623;
        c2[1, 1] := 0.1021101032532951;
        c2[1, 2] := 0.6365758882240111;
        c2[2, 1] := 0.1013729819392774;
        c2[2, 2] := 0.6267495975736321;
        c2[3, 1] := 0.1001476175660628;
        c2[3, 2] := 0.6104876178266819;
        c2[4, 1] := 0.9843854640428316e-1;
        c2[4, 2] := 0.5879603139195113;
        c2[5, 1] := 0.9625164534591696e-1;
        c2[5, 2] := 0.5594012291050210;
        c2[6, 1] := 0.9359356960417668e-1;
        c2[6, 2] := 0.5251016150410664;
        c2[7, 1] := 0.9047086748649986e-1;
        c2[7, 2] := 0.4854024475590397;
        c2[8, 1] := 0.8688856407189167e-1;
        c2[8, 2] := 0.4406826457109709;
        c2[9, 1] := 0.8284779224069856e-1;
        c2[9, 2] := 0.3913408089298914;
        c2[10, 1] := 0.7834154620997181e-1;
        c2[10, 2] := 0.3377643999400627;
        c2[11, 1] := 0.7334628941928766e-1;
        c2[11, 2] := 0.2802710651919946;
        c2[12, 1] := 0.6780290487362146e-1;
        c2[12, 2] := 0.2189770008083379;
        c2[13, 1] := 0.6156321231528423e-1;
        c2[13, 2] := 0.1534235999306070;
        c2[14, 1] := 0.5416797446761512e-1;
        c2[14, 2] := 0.8098664736760292e-1;
      elseif order == 30 then
        alpha := 0.1571200296252450;
        c2[1, 1] := 0.9908074847842124e-1;
        c2[1, 2] := 0.6289618807831557;
        c2[2, 1] := 0.9863509708328196e-1;
        c2[2, 2] := 0.6229164525571278;
        c2[3, 1] := 0.9774542692037148e-1;
        c2[3, 2] := 0.6108853364240036;
        c2[4, 1] := 0.9641490581986484e-1;
        c2[4, 2] := 0.5929869253412513;
        c2[5, 1] := 0.9464802912225441e-1;
        c2[5, 2] := 0.5693960175547550;
        c2[6, 1] := 0.9245027206218041e-1;
        c2[6, 2] := 0.5403402396359503;
        c2[7, 1] := 0.8982754584112941e-1;
        c2[7, 2] := 0.5060948065875106;
        c2[8, 1] := 0.8678535291732599e-1;
        c2[8, 2] := 0.4669749797983789;
        c2[9, 1] := 0.8332744242052199e-1;
        c2[9, 2] := 0.4233249626334694;
        c2[10, 1] := 0.7945356393775309e-1;
        c2[10, 2] := 0.3755006094498054;
        c2[11, 1] := 0.7515543969833788e-1;
        c2[11, 2] := 0.3238400339292700;
        c2[12, 1] := 0.7040879901685638e-1;
        c2[12, 2] := 0.2686072427439079;
        c2[13, 1] := 0.6515528854010540e-1;
        c2[13, 2] := 0.2098650589782619;
        c2[14, 1] := 0.5925168237177876e-1;
        c2[14, 2] := 0.1471138832654873;
        c2[15, 1] := 0.5225913954211672e-1;
        c2[15, 2] := 0.7775248839507864e-1;
      elseif order == 31 then
        alpha := 0.1545067022920929;
        c1[1] := 0.3100206996451866;
        c2[1, 1] := 0.9591020358831668e-1;
        c2[1, 2] := 0.6172474793293396;
        c2[2, 1] := 0.9530301275601203e-1;
        c2[2, 2] := 0.6088916323460413;
        c2[3, 1] := 0.9429332655402368e-1;
        c2[3, 2] := 0.5950511595503025;
        c2[4, 1] := 0.9288445429894548e-1;
        c2[4, 2] := 0.5758534119053522;
        c2[5, 1] := 0.9108073420087422e-1;
        c2[5, 2] := 0.5514734636081183;
        c2[6, 1] := 0.8888719137536870e-1;
        c2[6, 2] := 0.5221306199481831;
        c2[7, 1] := 0.8630901440239650e-1;
        c2[7, 2] := 0.4880834248148061;
        c2[8, 1] := 0.8335074993373294e-1;
        c2[8, 2] := 0.4496225358496770;
        c2[9, 1] := 0.8001502494376102e-1;
        c2[9, 2] := 0.4070602306679052;
        c2[10, 1] := 0.7630041338037624e-1;
        c2[10, 2] := 0.3607139804818122;
        c2[11, 1] := 0.7219760885744920e-1;
        c2[11, 2] := 0.3108783301229550;
        c2[12, 1] := 0.6768185077153345e-1;
        c2[12, 2] := 0.2577706252514497;
        c2[13, 1] := 0.6269571766328638e-1;
        c2[13, 2] := 0.2014081375889921;
        c2[14, 1] := 0.5710081766945065e-1;
        c2[14, 2] := 0.1412581515841926;
        c2[15, 1] := 0.5047740914807019e-1;
        c2[15, 2] := 0.7474725873250158e-1;
      elseif order == 32 then
        alpha := 0.1520196210848210;
        c2[1, 1] := 0.9322163554339406e-1;
        c2[1, 2] := 0.6101488690506050;
        c2[2, 1] := 0.9285233997694042e-1;
        c2[2, 2] := 0.6049832320721264;
        c2[3, 1] := 0.9211494244473163e-1;
        c2[3, 2] := 0.5946969295569034;
        c2[4, 1] := 0.9101176786042449e-1;
        c2[4, 2] := 0.5793791854364477;
        c2[5, 1] := 0.8954614071360517e-1;
        c2[5, 2] := 0.5591619969234026;
        c2[6, 1] := 0.8772216763680164e-1;
        c2[6, 2] := 0.5342177994699602;
        c2[7, 1] := 0.8554440426912734e-1;
        c2[7, 2] := 0.5047560942986598;
        c2[8, 1] := 0.8301735302045588e-1;
        c2[8, 2] := 0.4710187048140929;
        c2[9, 1] := 0.8014469519188161e-1;
        c2[9, 2] := 0.4332730387207936;
        c2[10, 1] := 0.7692807528893225e-1;
        c2[10, 2] := 0.3918021436411035;
        c2[11, 1] := 0.7336507157284898e-1;
        c2[11, 2] := 0.3468890521471250;
        c2[12, 1] := 0.6944555312763458e-1;
        c2[12, 2] := 0.2987898029050460;
        c2[13, 1] := 0.6514446669420571e-1;
        c2[13, 2] := 0.2476810747407199;
        c2[14, 1] := 0.6040544477732702e-1;
        c2[14, 2] := 0.1935412053397663;
        c2[15, 1] := 0.5509478650672775e-1;
        c2[15, 2] := 0.1358108994174911;
        c2[16, 1] := 0.4881064725720192e-1;
        c2[16, 2] := 0.7194819894416505e-1;
      elseif order == 33 then
        alpha := 0.1496489351138032;
        c1[1] := 0.3009752799176432;
        c2[1, 1] := 0.9041725460994505e-1;
        c2[1, 2] := 0.5995521047364046;
        c2[2, 1] := 0.8991117804113002e-1;
        c2[2, 2] := 0.5923764112099496;
        c2[3, 1] := 0.8906941547422532e-1;
        c2[3, 2] := 0.5804822013853129;
        c2[4, 1] := 0.8789442491445575e-1;
        c2[4, 2] := 0.5639663528946501;
        c2[5, 1] := 0.8638945831033775e-1;
        c2[5, 2] := 0.5429623519607796;
        c2[6, 1] := 0.8455834602616358e-1;
        c2[6, 2] := 0.5176379938389326;
        c2[7, 1] := 0.8240517431382334e-1;
        c2[7, 2] := 0.4881921474066189;
        c2[8, 1] := 0.7993380417355076e-1;
        c2[8, 2] := 0.4548502528082586;
        c2[9, 1] := 0.7714713890732801e-1;
        c2[9, 2] := 0.4178579388038483;
        c2[10, 1] := 0.7404596598181127e-1;
        c2[10, 2] := 0.3774715722484659;
        c2[11, 1] := 0.7062702339160462e-1;
        c2[11, 2] := 0.3339432938810453;
        c2[12, 1] := 0.6687952672391507e-1;
        c2[12, 2] := 0.2874950693388235;
        c2[13, 1] := 0.6277828912909767e-1;
        c2[13, 2] := 0.2382680702894708;
        c2[14, 1] := 0.5826808305383988e-1;
        c2[14, 2] := 0.1862073169968455;
        c2[15, 1] := 0.5321974125363517e-1;
        c2[15, 2] := 0.1307323751236313;
        c2[16, 1] := 0.4724820282032780e-1;
        c2[16, 2] := 0.6933542082177094e-1;
      elseif order == 34 then
        alpha := 0.1473858373968463;
        c2[1, 1] := 0.8801537152275983e-1;
        c2[1, 2] := 0.5929204288972172;
        c2[2, 1] := 0.8770594341007476e-1;
        c2[2, 2] := 0.5884653382247518;
        c2[3, 1] := 0.8708797598072095e-1;
        c2[3, 2] := 0.5795895850253119;
        c2[4, 1] := 0.8616320590689187e-1;
        c2[4, 2] := 0.5663615383647170;
        c2[5, 1] := 0.8493413175570858e-1;
        c2[5, 2] := 0.5488825092350877;
        c2[6, 1] := 0.8340387368687513e-1;
        c2[6, 2] := 0.5272851839324592;
        c2[7, 1] := 0.8157596213131521e-1;
        c2[7, 2] := 0.5017313864372913;
        c2[8, 1] := 0.7945402670834270e-1;
        c2[8, 2] := 0.4724089864574216;
        c2[9, 1] := 0.7704133559556429e-1;
        c2[9, 2] := 0.4395276256463053;
        c2[10, 1] := 0.7434009635219704e-1;
        c2[10, 2] := 0.4033126590648964;
        c2[11, 1] := 0.7135035113853376e-1;
        c2[11, 2] := 0.3639961488919042;
        c2[12, 1] := 0.6806813160738834e-1;
        c2[12, 2] := 0.3218025212900124;
        c2[13, 1] := 0.6448214312000864e-1;
        c2[13, 2] := 0.2769235521088158;
        c2[14, 1] := 0.6056719318430530e-1;
        c2[14, 2] := 0.2294693573271038;
        c2[15, 1] := 0.5626925196925040e-1;
        c2[15, 2] := 0.1793564218840015;
        c2[16, 1] := 0.5146352031547277e-1;
        c2[16, 2] := 0.1259877129326412;
        c2[17, 1] := 0.4578069074410591e-1;
        c2[17, 2] := 0.6689147319568768e-1;
      elseif order == 35 then
        alpha := 0.1452224267615486;
        c1[1] := 0.2926764667564367;
        c2[1, 1] := 0.8551731299267280e-1;
        c2[1, 2] := 0.5832758214629523;
        c2[2, 1] := 0.8509109732853060e-1;
        c2[2, 2] := 0.5770596582643844;
        c2[3, 1] := 0.8438201446671953e-1;
        c2[3, 2] := 0.5667497616665494;
        c2[4, 1] := 0.8339191981579831e-1;
        c2[4, 2] := 0.5524209816238369;
        c2[5, 1] := 0.8212328610083385e-1;
        c2[5, 2] := 0.5341766459916322;
        c2[6, 1] := 0.8057906332198853e-1;
        c2[6, 2] := 0.5121470053512750;
        c2[7, 1] := 0.7876247299954955e-1;
        c2[7, 2] := 0.4864870722254752;
        c2[8, 1] := 0.7667670879950268e-1;
        c2[8, 2] := 0.4573736721705665;
        c2[9, 1] := 0.7432449556218945e-1;
        c2[9, 2] := 0.4250013835198991;
        c2[10, 1] := 0.7170742126011575e-1;
        c2[10, 2] := 0.3895767735915445;
        c2[11, 1] := 0.6882488171701314e-1;
        c2[11, 2] := 0.3513097926737368;
        c2[12, 1] := 0.6567231746957568e-1;
        c2[12, 2] := 0.3103999917596611;
        c2[13, 1] := 0.6223804362223595e-1;
        c2[13, 2] := 0.2670123611280899;
        c2[14, 1] := 0.5849696460782910e-1;
        c2[14, 2] := 0.2212298104867592;
        c2[15, 1] := 0.5439628409499822e-1;
        c2[15, 2] := 0.1729443731341637;
        c2[16, 1] := 0.4981540179136920e-1;
        c2[16, 2] := 0.1215462157134930;
        c2[17, 1] := 0.4439981033536435e-1;
        c2[17, 2] := 0.6460098363520967e-1;
      elseif order == 36 then
        alpha := 0.1431515914458580;
        c2[1, 1] := 0.8335881847130301e-1;
        c2[1, 2] := 0.5770670512160201;
        c2[2, 1] := 0.8309698922852212e-1;
        c2[2, 2] := 0.5731929100172432;
        c2[3, 1] := 0.8257400347039723e-1;
        c2[3, 2] := 0.5654713811993058;
        c2[4, 1] := 0.8179117911600136e-1;
        c2[4, 2] := 0.5539556343603020;
        c2[5, 1] := 0.8075042173126963e-1;
        c2[5, 2] := 0.5387245649546684;
        c2[6, 1] := 0.7945413151258206e-1;
        c2[6, 2] := 0.5198817177723069;
        c2[7, 1] := 0.7790506514288866e-1;
        c2[7, 2] := 0.4975537629595409;
        c2[8, 1] := 0.7610613635339480e-1;
        c2[8, 2] := 0.4718884193866789;
        c2[9, 1] := 0.7406012816626425e-1;
        c2[9, 2] := 0.4430516443136726;
        c2[10, 1] := 0.7176927060205631e-1;
        c2[10, 2] := 0.4112237708115829;
        c2[11, 1] := 0.6923460172504251e-1;
        c2[11, 2] := 0.3765940116389730;
        c2[12, 1] := 0.6645495833489556e-1;
        c2[12, 2] := 0.3393522147815403;
        c2[13, 1] := 0.6342528888937094e-1;
        c2[13, 2] := 0.2996755899575573;
        c2[14, 1] := 0.6013361864949449e-1;
        c2[14, 2] := 0.2577053294053830;
        c2[15, 1] := 0.5655503081322404e-1;
        c2[15, 2] := 0.2135004731531631;
        c2[16, 1] := 0.5263798119559069e-1;
        c2[16, 2] := 0.1669320999865636;
        c2[17, 1] := 0.4826589873626196e-1;
        c2[17, 2] := 0.1173807590715484;
        c2[18, 1] := 0.4309819397289806e-1;
        c2[18, 2] := 0.6245036108880222e-1;
      elseif order == 37 then
        alpha := 0.1411669104782917;
        c1[1] := 0.2850271036215707;
        c2[1, 1] := 0.8111958235023328e-1;
        c2[1, 2] := 0.5682412610563970;
        c2[2, 1] := 0.8075727567979578e-1;
        c2[2, 2] := 0.5628142923227016;
        c2[3, 1] := 0.8015440554413301e-1;
        c2[3, 2] := 0.5538087696879930;
        c2[4, 1] := 0.7931239302677386e-1;
        c2[4, 2] := 0.5412833323304460;
        c2[5, 1] := 0.7823314328639347e-1;
        c2[5, 2] := 0.5253190555393968;
        c2[6, 1] := 0.7691895211595101e-1;
        c2[6, 2] := 0.5060183741977191;
        c2[7, 1] := 0.7537237072011853e-1;
        c2[7, 2] := 0.4835036020049034;
        c2[8, 1] := 0.7359601294804538e-1;
        c2[8, 2] := 0.4579149413954837;
        c2[9, 1] := 0.7159227884849299e-1;
        c2[9, 2] := 0.4294078049978829;
        c2[10, 1] := 0.6936295002846032e-1;
        c2[10, 2] := 0.3981491350382047;
        c2[11, 1] := 0.6690857785828917e-1;
        c2[11, 2] := 0.3643121502867948;
        c2[12, 1] := 0.6422751692085542e-1;
        c2[12, 2] := 0.3280684291406284;
        c2[13, 1] := 0.6131430866206096e-1;
        c2[13, 2] := 0.2895750997170303;
        c2[14, 1] := 0.5815677249570920e-1;
        c2[14, 2] := 0.2489521814805720;
        c2[15, 1] := 0.5473023527947980e-1;
        c2[15, 2] := 0.2062377435955363;
        c2[16, 1] := 0.5098441033167034e-1;
        c2[16, 2] := 0.1612849131645336;
        c2[17, 1] := 0.4680658811093562e-1;
        c2[17, 2] := 0.1134672937045305;
        c2[18, 1] := 0.4186928031694695e-1;
        c2[18, 2] := 0.6042754777339966e-1;
      elseif order == 38 then
        alpha := 0.1392625697140030;
        c2[1, 1] := 0.7916943373658329e-1;
        c2[1, 2] := 0.5624158631591745;
        c2[2, 1] := 0.7894592250257840e-1;
        c2[2, 2] := 0.5590219398777304;
        c2[3, 1] := 0.7849941672384930e-1;
        c2[3, 2] := 0.5522551628416841;
        c2[4, 1] := 0.7783093084875645e-1;
        c2[4, 2] := 0.5421574325808380;
        c2[5, 1] := 0.7694193770482690e-1;
        c2[5, 2] := 0.5287909941093643;
        c2[6, 1] := 0.7583430534712885e-1;
        c2[6, 2] := 0.5122376814029880;
        c2[7, 1] := 0.7451020436122948e-1;
        c2[7, 2] := 0.4925978555548549;
        c2[8, 1] := 0.7297197617673508e-1;
        c2[8, 2] := 0.4699889739625235;
        c2[9, 1] := 0.7122194706992953e-1;
        c2[9, 2] := 0.4445436860615774;
        c2[10, 1] := 0.6926216260386816e-1;
        c2[10, 2] := 0.4164072786327193;
        c2[11, 1] := 0.6709399961255503e-1;
        c2[11, 2] := 0.3857341621868851;
        c2[12, 1] := 0.6471757977022456e-1;
        c2[12, 2] := 0.3526828388476838;
        c2[13, 1] := 0.6213084287116965e-1;
        c2[13, 2] := 0.3174082831364342;
        c2[14, 1] := 0.5932799638550641e-1;
        c2[14, 2] := 0.2800495563550299;
        c2[15, 1] := 0.5629672408524944e-1;
        c2[15, 2] := 0.2407078154782509;
        c2[16, 1] := 0.5301264751544952e-1;
        c2[16, 2] := 0.1994026830553859;
        c2[17, 1] := 0.4942673259817896e-1;
        c2[17, 2] := 0.1559719194038917;
        c2[18, 1] := 0.4542996716979947e-1;
        c2[18, 2] := 0.1097844277878470;
        c2[19, 1] := 0.4070720755433961e-1;
        c2[19, 2] := 0.5852181110523043e-1;
      elseif order == 39 then
        alpha := 0.1374332900196804;
        c1[1] := 0.2779468246419593;
        c2[1, 1] := 0.7715084161825772e-1;
        c2[1, 2] := 0.5543001331300056;
        c2[2, 1] := 0.7684028301163326e-1;
        c2[2, 2] := 0.5495289890712267;
        c2[3, 1] := 0.7632343924866024e-1;
        c2[3, 2] := 0.5416083298429741;
        c2[4, 1] := 0.7560141319808483e-1;
        c2[4, 2] := 0.5305846713929198;
        c2[5, 1] := 0.7467569064745969e-1;
        c2[5, 2] := 0.5165224112570647;
        c2[6, 1] := 0.7354807648551346e-1;
        c2[6, 2] := 0.4995030679271456;
        c2[7, 1] := 0.7222060351121389e-1;
        c2[7, 2] := 0.4796242430956156;
        c2[8, 1] := 0.7069540462458585e-1;
        c2[8, 2] := 0.4569982440368368;
        c2[9, 1] := 0.6897453353492381e-1;
        c2[9, 2] := 0.4317502624832354;
        c2[10, 1] := 0.6705970959388781e-1;
        c2[10, 2] := 0.4040159353969854;
        c2[11, 1] := 0.6495194541066725e-1;
        c2[11, 2] := 0.3739379843169939;
        c2[12, 1] := 0.6265098412417610e-1;
        c2[12, 2] := 0.3416613843816217;
        c2[13, 1] := 0.6015440984955930e-1;
        c2[13, 2] := 0.3073260166338746;
        c2[14, 1] := 0.5745615876877304e-1;
        c2[14, 2] := 0.2710546723961181;
        c2[15, 1] := 0.5454383762391338e-1;
        c2[15, 2] := 0.2329316824061170;
        c2[16, 1] := 0.5139340231935751e-1;
        c2[16, 2] := 0.1929604256043231;
        c2[17, 1] := 0.4795705862458131e-1;
        c2[17, 2] := 0.1509655259246037;
        c2[18, 1] := 0.4412933231935506e-1;
        c2[18, 2] := 0.1063130748962878;
        c2[19, 1] := 0.3960672309405603e-1;
        c2[19, 2] := 0.5672356837211527e-1;
      elseif order == 40 then
        alpha := 0.1356742655825434;
        c2[1, 1] := 0.7538038374294594e-1;
        c2[1, 2] := 0.5488228264329617;
        c2[2, 1] := 0.7518806529402738e-1;
        c2[2, 2] := 0.5458297722483311;
        c2[3, 1] := 0.7480383050347119e-1;
        c2[3, 2] := 0.5398604576730540;
        c2[4, 1] := 0.7422847031965465e-1;
        c2[4, 2] := 0.5309482987446206;
        c2[5, 1] := 0.7346313704205006e-1;
        c2[5, 2] := 0.5191429845322307;
        c2[6, 1] := 0.7250930053201402e-1;
        c2[6, 2] := 0.5045099368431007;
        c2[7, 1] := 0.7136868456879621e-1;
        c2[7, 2] := 0.4871295553902607;
        c2[8, 1] := 0.7004317764946634e-1;
        c2[8, 2] := 0.4670962098860498;
        c2[9, 1] := 0.6853470921527828e-1;
        c2[9, 2] := 0.4445169164956202;
        c2[10, 1] := 0.6684507689945471e-1;
        c2[10, 2] := 0.4195095960479698;
        c2[11, 1] := 0.6497570123412630e-1;
        c2[11, 2] := 0.3922007419030645;
        c2[12, 1] := 0.6292726794917847e-1;
        c2[12, 2] := 0.3627221993494397;
        c2[13, 1] := 0.6069918741663154e-1;
        c2[13, 2] := 0.3312065181294388;
        c2[14, 1] := 0.5828873983769410e-1;
        c2[14, 2] := 0.2977798532686911;
        c2[15, 1] := 0.5568964389813015e-1;
        c2[15, 2] := 0.2625503293999835;
        c2[16, 1] := 0.5288947816690705e-1;
        c2[16, 2] := 0.2255872486520188;
        c2[17, 1] := 0.4986456327645859e-1;
        c2[17, 2] := 0.1868796731919594;
        c2[18, 1] := 0.4656832613054458e-1;
        c2[18, 2] := 0.1462410193532463;
        c2[19, 1] := 0.4289867647614935e-1;
        c2[19, 2] := 0.1030361558710747;
        c2[20, 1] := 0.3856310684054106e-1;
        c2[20, 2] := 0.5502423832293889e-1;
      elseif order == 41 then
        alpha := 0.1339811106984253;
        c1[1] := 0.2713685065531391;
        c2[1, 1] := 0.7355140275160984e-1;
        c2[1, 2] := 0.5413274778282860;
        c2[2, 1] := 0.7328319082267173e-1;
        c2[2, 2] := 0.5371064088294270;
        c2[3, 1] := 0.7283676160772547e-1;
        c2[3, 2] := 0.5300963437270770;
        c2[4, 1] := 0.7221298133014343e-1;
        c2[4, 2] := 0.5203345998371490;
        c2[5, 1] := 0.7141302173623395e-1;
        c2[5, 2] := 0.5078728971879841;
        c2[6, 1] := 0.7043831559982149e-1;
        c2[6, 2] := 0.4927768111819803;
        c2[7, 1] := 0.6929049381827268e-1;
        c2[7, 2] := 0.4751250308594139;
        c2[8, 1] := 0.6797129849758392e-1;
        c2[8, 2] := 0.4550083840638406;
        c2[9, 1] := 0.6648246325101609e-1;
        c2[9, 2] := 0.4325285673076087;
        c2[10, 1] := 0.6482554675958526e-1;
        c2[10, 2] := 0.4077964789091151;
        c2[11, 1] := 0.6300169683004558e-1;
        c2[11, 2] := 0.3809299858742483;
        c2[12, 1] := 0.6101130648543355e-1;
        c2[12, 2] := 0.3520508315700898;
        c2[13, 1] := 0.5885349417435808e-1;
        c2[13, 2] := 0.3212801560701271;
        c2[14, 1] := 0.5652528148656809e-1;
        c2[14, 2] := 0.2887316252774887;
        c2[15, 1] := 0.5402021575818373e-1;
        c2[15, 2] := 0.2545001287790888;
        c2[16, 1] := 0.5132588802608274e-1;
        c2[16, 2] := 0.2186415296842951;
        c2[17, 1] := 0.4841900639702602e-1;
        c2[17, 2] := 0.1811322622296060;
        c2[18, 1] := 0.4525419574485134e-1;
        c2[18, 2] := 0.1417762065404688;
        c2[19, 1] := 0.4173260173087802e-1;
        c2[19, 2] := 0.9993834530966510e-1;
        c2[20, 1] := 0.3757210572966463e-1;
        c2[20, 2] := 0.5341611499960143e-1;
      else
        Streams.error("Input argument order (= " + String(order) +
          ") of Bessel filter is not in the range 1..41");
      end if;

      annotation (Documentation(info="<html>
<p>
The transfer function H(p) of a <i>n</i> 'th order Bessel filter is given by
</p>
<blockquote><pre>
        Bn(0)
H(p) = -------
        Bn(p)
</pre></blockquote>
<p>
with the denominator polynomial
</p>
<blockquote><pre>
         n             n  (2n - k)!       p^k
Bn(p) = sum c_k*p^k = sum ----------- * -------   (1)
        k=0           k=0 (n - k)!k!    2^(n-k)
</pre></blockquote>
<p>
and the numerator
</p>
<blockquote><pre>
               (2n)!     1
Bn(0) = c_0 = ------- * ---- .                     (2)
                n!      2^n
</pre></blockquote>
<p>
Although the coefficients c_k are integer numbers, it is not advisable to use the
polynomials in an unfactorized form because the coefficients are fast growing with order
n (c_0 is approximately 0.3e24 and 0.8e59 for order n=20 and order n=40
respectively).
</p>
<p>
Therefore, the polynomial Bn(p) is factorized to first and second order polynomials with
real coefficients corresponding to zeros and poles representation that is used in this library.
</p>
<p>
The function returns the coefficients which resulted from factorization of the normalized transfer function
</p>
<blockquote><pre>
H'(p') = H(p),  p' = p/w0
</pre></blockquote>
<p>
as well as
</p>
<blockquote><pre>
alpha = 1/w0
</pre></blockquote>
<p>
the reciprocal of the cut of frequency w0 where the gain of the transfer function is
decreased 3dB.
</p>
<p>
Both, coefficients and cut off frequency were calculated symbolically and were eventually evaluated
with high precision calculation. The results were stored in this function as real
numbers.
</p>

<h4>Calculation of normalized Bessel filter coefficients</h4>
<p>
Equation
<blockquote><pre>
abs(H(j*w0)) = abs(Bn(0)/Bn(j*w0)) = 10^(-3/20)
</pre></blockquote>
<p>
which must be fulfilled for cut off frequency w = w0 leads to
</p>
<blockquote><pre>
[Re(Bn(j*w0))]^2 + [Im(Bn(j*w0))]^2 - (Bn(0)^2)*10^(3/10) = 0
</pre></blockquote>
<p>
which has exactly one real solution w0 for each order n. This solutions of w0 are
calculated symbolically first and evaluated by using high precise values of the
coefficients c_k calculated by following (1) and (2).
</p>
<p>
With w0, the coefficients of the factorized polynomial can be computed by calculating the
zeros of the denominator polynomial
</p>
<blockquote><pre>
        n
Bn(p) = sum w0^k*c_k*(p/w0)^k
        k=0
</pre></blockquote>
<p>
of the normalized transfer function H'(p'). There exist n/2 of conjugate complex
pairs of zeros (beta +-j*gamma) if n is even and one additional real zero (alpha) if n is
odd. Finally, the coefficients a, b1_k, b2_k of the polynomials
</p>
<blockquote><pre>
a*p + 1,  n is odd
</pre></blockquote>
<p>
and
</p>
<blockquote><pre>
b2_k*p^2 + b1_k*p + 1,   k = 1,... div(n,2)
</pre></blockquote>
<p>
results from
<blockquote><pre>
a = -1/alpha
</pre></blockquote>
<p>
and
</p>
<blockquote><pre>
b2_k = 1/(beta_k^2 + gamma_k^2) b1_k = -2*beta_k/(beta_k^2 + gamma_k^2)
</pre></blockquote>
</html>"));
    end BesselCoefficients;

    function checkRepresentation
      "Check whether the system on file is represented by zeros and poles (z, p) or first and second order polynomials (n1, n2, d1, d2)"
      import Modelica_LinearSystems2.ZerosAndPoles;
      input String fileName="zp.mat" "Name of the zeros and poles data file"
        annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
          caption="state space system data file")));
      output Boolean iszp=true;
    protected
      Integer m=0;

    algorithm
      m := ZerosAndPoles.Internal.findMatrixName(fileName, "z");
      m := m + ZerosAndPoles.Internal.findMatrixName(fileName, "p");
      iszp := m == 2;

      annotation (Documentation(info="<html>
<p>
The function output is true if the system is given in zeros and poles representation.
Therefore, it is assumend that the used array names are &quot;z&quot; and &quot;p&quot;
or &quot;n1, n2, d1&quot; and &quot;d2&quot; respectively.
</p>
</html>"));
    end checkRepresentation;

  encapsulated function filter
      "Generate the data record of a ZerosAndPoles transfer function from a filter description"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Utilities.Types;
      import Modelica_LinearSystems2.ZerosAndPoles;

      import MMath = Modelica.Math;

      input Modelica_LinearSystems2.Utilities.Types.AnalogFilter analogFilter=Modelica_LinearSystems2.Utilities.Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
      input Modelica_LinearSystems2.Utilities.Types.FilterType filterType=Modelica_LinearSystems2.Utilities.Types.FilterType.LowPass "Type of filter (LowPass/HighPass)";
    input Integer order(min=1) = 2 "Order of filter";
    input Modelica.SIunits.Frequency f_cut=1/(2*Modelica.Constants.pi)
        "Cut-off frequency (default is w_cut = 1 rad/s)";
    input Real gain=1.0
        "Gain (= amplitude of frequency response at zero frequency)";
    input Real A_ripple(unit="dB") = 0.5
        "Pass band ripple for Chebyshev filter (otherwise not used)";
    input Boolean normalized=true
        "True, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
    output ZerosAndPoles.Internal.ZerosAndPoles filter(
      redeclare Real n1[if filterType == Types.FilterType.LowPass then 0 else (
        if analogFilter == Types.AnalogFilter.CriticalDamping then order else
        mod(order, 2))],
      redeclare Real n2[if filterType == Types.FilterType.LowPass then 0 else (
        if analogFilter == Types.AnalogFilter.CriticalDamping then 0 else
        integer(order/2)),2],
      redeclare Real d1[if analogFilter == Types.AnalogFilter.CriticalDamping then
              order else mod(order, 2)],
      redeclare Real d2[if analogFilter == Types.AnalogFilter.CriticalDamping then
              0 else integer(order/2),2]) "Filter transfer function";

    protected
    Integer n_num1=size(filter.n1, 1);
    Integer n_num2=size(filter.n2, 1);
    Integer n_den1=size(filter.d1, 1);
    Integer n_den2=size(filter.d2, 1);
    Integer n_num=n_num1 + 2*n_num2;
    Integer n_den=n_den1 + 2*n_den2;
    Real pi=Modelica.Constants.pi;
    Boolean evenOrder=mod(order, 2) == 0
        "True, if even filter order, otherwise uneven";
    Modelica.SIunits.AngularVelocity w_cut=2*pi*f_cut
        "Cut-off angular frequency";
    Real w_cut2 "= w_cut*w_cut";
    Real alpha=1.0 "Frequency correction factor";
    Real alpha2 "= alpha*alpha";

    Real alphax;

    Real epsilon "Ripple size";
    Real fac "arsinh(epsilon)";
    Real A2 "poleReal^2 + poleImag^2";
    Real A "Amplitude at w_cut";
    Real num1[n_num1]
        "[p] coefficients of numerator first order polynomials (a*p + 1)";
    Real num2[n_num2,2]
        "[p^2, p] coefficients of numerator second order polynomials (b*p^2 + a*p + 1I)";
    Real den1[n_den1]
        "[p] coefficients of denominator first order polynomials (a*p + 1)";
    Real den2[n_den2,2]
        "[p^2, p] coefficients of denominator second order polynomials (b*p^2 + a*p + 1I)";
    Real aux;
    Real k;
  algorithm
    // Set properties that are common for all filters
    filter.k := gain;

    /* Compute filter coefficients of prototype low pass filter. If another filter
     characteristics is desired (e.g. high pass filter), it is derived
     from the low pass filter coefficients below
  */
    if analogFilter == Types.AnalogFilter.CriticalDamping then
      if normalized then
        alpha := sqrt(2^(1/order) - 1);
  //alpha := sqrt(10^(3/10/order)-1)
      else
        alpha := 1;
      end if;
      for i in 1:n_den1 loop
        den1[i] := alpha;
      end for;

    elseif analogFilter == Types.AnalogFilter.Bessel then
      (den1,den2,alpha) := ZerosAndPoles.Internal.BesselCoefficients(order);
      if not normalized then
        alpha2 := alpha*alpha;
        for i in 1:n_den2 loop
          den2[i, 1] := den2[i, 1]*alpha2;
          den2[i, 2] := den2[i, 2]*alpha;
        end for;
        if not evenOrder then
          den1[1] := den1[1]*alpha;
        end if;
      end if;

    elseif analogFilter == Types.AnalogFilter.Butterworth then
       // Original filter is already normalized
      for i in 1:n_den2 loop
        den2[i, 1] := 1.0;
        den2[i, 2] := -2*cos(pi*(0.5 + (i - 0.5)/order));
      end for;
      if not evenOrder then
        den1[1] := 1.0;
      end if;

    elseif analogFilter == Types.AnalogFilter.Chebyshev then
      epsilon := sqrt(10^(A_ripple/10) - 1);
      fac := MMath.asinh(1/epsilon)/order;

      if evenOrder then
         for i in 1:n_den2 loop
            den2[i,1] :=1/(cosh(fac)^2 - cos((2*i - 1)*pi/(2*order))^2);
            den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos((2*i - 1)*pi/(2*order));
         end for;
      else
         den1[1] := 1/sinh(fac);
         for i in 1:n_den2 loop
            den2[i,1] :=1/(cosh(fac)^2 - cos(i*pi/order)^2);
            den2[i,2] :=2*den2[i, 1]*sinh(fac)*cos(i*pi/order);
         end for;
      end if;

       /* Transformation of filter transfer function with "new(p) = alpha*p"
        in order that the filter transfer function has an amplitude of
        1/sqrt(2) at the cutoff frequency
     */
      if normalized then
        alpha := ZerosAndPoles.Internal.normalizationFactor(den1, den2);
        alpha2 := alpha*alpha;
        for i in 1:n_den2 loop
          den2[i, 1] := den2[i, 1]*alpha2;
          den2[i, 2] := den2[i, 2]*alpha;
        end for;
        if not evenOrder then
          den1[1] := den1[1]*alpha;
        end if;
      end if;

    else
      Streams.error("analogFilter (= " + String(analogFilter) +
        ") is not supported");
    end if;

    // Compute amplitude at w=1
  /*
  A := 1.0;
  for i in 1:n_den2 loop
     A := A*(1 + den2[i,2]^2 - 2*den2[i,1]
               + den2[i,1]^2);
  end for;
  for i in 1:n_den1 loop
     A := A*(1 + den1[i]^2);
  end for;
  A := 1/sqrt(A);
  Streams.print("A = " + String(A));
*/

    // Determine normalized denominator polynomials with highest power of p equal to one
    filter.n1 := zeros(n_num1);
    filter.n2 := zeros(n_num2, 2);
    (filter.d1,filter.d2,k) := ZerosAndPoles.Internal.filterToNormalized(den1, den2);
    filter.k := filter.k/k;

    // Compute desired filter characteristics from low pass filter coefficients
    if filterType == Types.FilterType.HighPass then
       /* The high pass filter is derived from the low pass filter by
        the transformation new(p) = 1/p
        1/(p^2 + a*p + b) -> 1/((1/p)^2 + a*(1/p) + b) = (1/b)*p^2 / (p^2 + (a/b)*p + 1/b)
        1/(p + a)         -> 1/((1/p) + a) = (1/a)*p / (p + (1/a))
     */
      assert(n_num1 == n_den1 and n_num2 == n_den2,
        "Internal error 1, should not occur");
      filter.n1 := zeros(n_num1);
      filter.n2 := zeros(n_num2, 2);
      for i in 1:n_num1 loop
        filter.k := filter.k/filter.d1[i];
        filter.d1[i] := 1/filter.d1[i];
      end for;
      for i in 1:n_num2 loop
        filter.k := filter.k/filter.d2[i, 2];
        filter.d2[i, 1] := filter.d2[i, 1]/filter.d2[i, 2];
        filter.d2[i, 2] := 1/filter.d2[i, 2];
      end for;
    end if;

    /* Change filter coefficients according to transformation new(p) = p/w_cut
     Numerator  :     (p/w)^2 + a*(p/w) + b = (1/w^2)*(p^2 + (a*w)*p + b*w^2)
                                  (p/w) + a = (1/w)*(p + w*a)
     Denominator: 1/((p/w)^2 + a*(p/w) + b) = w^2/(p^2 + (a*w)*p + w^2/b)
                              1/((p/w) + a) = w/(p + w*a)
  */
    w_cut2 := w_cut*w_cut;
    filter.k := filter.k*w_cut^(n_den1 + 2*n_den2 - n_num1 - 2*n_num2);
    filter.n1 := w_cut*filter.n1;
    filter.d1 := w_cut*filter.d1;
    filter.n2 := [w_cut*filter.n2[:, 1],w_cut2*filter.n2[:, 2]];
    filter.d2 := [w_cut*filter.d2[:, 1],w_cut2*filter.d2[:, 2]];

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
filterFunction = ZerosAndPoles.Internal<b>filter</b>(
  analogFilter,
  filterType,
  order,
  f_cut,
  gain,
  A_ripple,
  normalized)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a ZerosAndPoles transfer function
description of low and high pass filters. For more details see also
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Literature\">[Tietze2002]</a>, pp. 815-852.
</p>
<p>
Typical frequency responses for the four supported low pass filter types
are shown in the next figure (this figure was generated with function
<a href=\"modelica://Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter2\">Examples.ZerosAndPoles.plotBodeFilter2</a>):
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4Filters.png\">
</p>
<p>
The step responses of the same low pass filters are shown in the next figure,
starting from a steady state initial filter with initial input = 0.2:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4FiltersStepResponse.png\">
</p>
<p>
Obviously, the frequency responses give a somewhat wrong impression
of the filter characteristics: Although Butterworth and Chebyshev
filters have a significantly steeper magnitude as the
CriticalDamping and Bessel filters, the step responses of
the latter ones are much better. This means for example, that
a CriticalDamping or a Bessel filter should be selected,
if a filter is mainly used to make a non-linear inverse model
realizable.
</p>

<p>
Typical frequency responses for the four supported high pass filter types
are shown in the next figure (generated with function
<a href=\"modelica://Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter3\">Examples.ZerosAndPoles.plotBodeFilter3</a>):
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4Filters.png\">
</p>
<p>
The corresponding step responses of these high pass filters are
shown in the next figure:
</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4FiltersStepResponse.png\">
</p>
<p>
All filters are available in <b>normalized</b> (default) and non-normalized form.
In the normalized form, the amplitude of the filter transfer function
at the cutoff frequency is 1/sqrt(2) (= 3 dB). Note, when comparing the filters
of this function with other software systems, the setting of &quot;normalized&quot;
has to be selected appropriately. For example, the signal processing
toolbox of Matlab provides the filters in non-normalized form and
therefore a comparison makes only sense, if normalized = <b>false</b>
is set.
</p>

<h4>Example</h4>
<blockquote><pre>
  Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping;
  Integer order=2;
  Modelica.SIunits.Frequency f_cut=10;

  ZerosAndPoles zp_filter;

<b>algorithm</b>
  zp_filter=Modelica_LinearSystems2.ZerosAndPoles.Design.filter(
    order=order,
    f_cut=f_cut,
    analogFilter=analogFilter);

// zp_filter = 9530.93/( (p + 97.6265)^2 )
</pre></blockquote>
</html>"));
  end filter;

    function filterToNormalized
      "Given [p^2,p] and [p] coefficients, transform to normalized form with highest power of p equal 1"

      input Real c1[:] "[p] coefficients of polynomials (c1[i]*p + 1)";
      input Real c2[:,2]
        "[p^2, p] coefficients of polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
      output Real n1[size(c1, 1)]
        "[p^0] coefficients of polynomials c1[i]*(p+1/c1[i])";
      output Real n2[size(c2, 1),2]
        "[p, p^0] coefficients of polynomials c2[i,1]*(p^2 + (c2[i,2]/c2[i,1])*p + (1/c2[i,1]))";
      output Real k "Gain (product(1/c1[i])*(1/c2[i,1])";
    algorithm
      k := 1.0;
      for i in 1:size(c1, 1) loop
        k := k*c1[i];
        n1[i] := 1/c1[i];
      end for;

      for i in 1:size(c2, 1) loop
        k := k*c2[i, 1];
        n2[i, 1] := c2[i, 2]/c2[i, 1];
        n2[i, 2] := 1/c2[i, 1];
      end for;
    end filterToNormalized;

    function findMatrixName
      "Find out whether matrix matName exists in file filename"
      input String filename;
      input String matName="z";

      output Integer m;

    external "C" m = findMatrixName(
            filename,
            matName,
            "NoClass");

      annotation (Include="#include <matrixop.h>
#include <matrixop.c>


#if !defined(DYMOLA_DSPACE) && !defined(NO_FILE)
#include <amat.h>
#include <sprwat.h>
#endif

extern int findMatrixName(const char* fil,const char* matname, char *noClass) {
int found=0;


#if !defined(DYMOLA_DSPACE) && !defined(NO_FILE)
{
        AmatGetFile afile;
        Amatrix amatrix;

        int ret=amatGetOpen((char*)fil,noClass,(char*)0,&afile);

        Assert(ret==0,amatError);
        for(;ret==0 && !found;)
{
                amatInit(&amatrix);
                ret=amatGetMatrix(&afile, &amatrix);
                if (ret<=1 && strcmp(matname,amatrix.name)==0)
                  found=1;
                else
                  found=0;

                amatDel(&amatrix);
        }
        amatGetClose(&afile);
}
#else
        Assert(false, 'nn');
#endif
        return found;
}");
    end findMatrixName;

    function firstOrderToString
      "Transform vector of coefficients of first order polynomials to a string representation"
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math;

      input Real c[:] = fill(0.0,0)
        "Coefficients of first order polynomials: polynom(p) = p + c[i]";
      input Integer significantDigits=6
        "Number of significant digits that are shown";
      input String name="p" "Independent variable name used for printing";
      input Boolean normalized=false
        "= true, the polynomials in the string are represented as p/c[i] + 1, provided c[i]<>0";
      output String s="";
      output Real gain=1.0
        "If normalized=true, the product(c[i]) (for i, where c[i]<>0), otherwise gain=1.0";
    protected
      Integer nc=size(c, 1);
      Real cs[nc];
      Real cc[nc];
      Integer i;
      Integer j;
      Integer j2;
      Integer nj;
      Integer k;
      constant Real smallNumber = 100*Modelica.Constants.small;
      constant Real eps = 10*Modelica.Constants.eps;
    algorithm
      // Change coefficients, if normalized output, and sort them
      for i in 1:nc loop
         if Math.isEqual(c[i], 0.0, smallNumber) then
            cs[i] := 0.0;
         elseif normalized then
            cs[i] := 1/c[i];
            gain := gain*c[i];
         else
            cs[i] := c[i];
         end if;
      end for;
      cc :=Modelica.Math.Vectors.sort(cs);

      // Move zeros to the beginning
      j :=0;
      k :=nc + 1;
      for i in nc:-1:1 loop
         if cc[i] == 0.0 then
            j := j + 1;
            cs[j] := 0.0;
         else
            k := k - 1;
            cs[k] := cc[i];
         end if;
      end for;

      // Transform coefficients to string
      if j == 1 then
         s := name;
      elseif j > 1 then
         s := name + "^" + String(j);
      end if;

      i :=j + 1;
      if normalized then
        while i <= nc loop
          if i > 1 then
            s := s + "*";
          end if;
          if Math.isEqual(cs[i], 1.0, eps) then
             s := s + "(";
          elseif Math.isEqual(cs[i], -1.0, eps) then
             s := s + "(-";
          else
             s := s + "(" + String(cs[i], significantDigits=significantDigits) + "*";
          end if;
          s :=s + name + " + 1)";
          j2 := sameVectorElements(cs, i);
          nj := j2 - i + 1;
          if nj > 1 then
            s := s + "^" + String(nj);
          end if;
          i := j2 + 1;
        end while;
      else
        while i <= nc loop
          if i > 1 then
            s := s + "*";
          end if;
          if cs[i] > 0 then
             s := s + "(" + name + " + ";
          else
             s := s + "(" + name + " - ";
          end if;
          s  := s + String(abs(cs[i]), significantDigits=significantDigits) + ")";
          j2 := sameVectorElements(cs, i);
          nj := j2 - i + 1;
          if nj > 1 then
            s := s + "^" + String(nj);
          end if;
          i := j2 + 1;
        end while;
      end if;
    end firstOrderToString;

    function frequencyRange "Determine min. and max. resonance frequencies"
      import Modelica;
      import Complex;

      input Real poly1[:];
      input Real poly2[:,2];
      output Boolean w_found=false;
      output Modelica.SIunits.AngularVelocity w_min;
      output Modelica.SIunits.AngularVelocity w_max;
    protected
      Integer order=size(poly1, 1) + 2*size(poly2, 1);
      Integer n_real=numberOfRealZeros(poly1, poly2);
      Real zeros1[n_real];
      Complex zeros2[:]=fill(Complex(0, 0), integer((order - n_real)/2));
      Real w;
    algorithm
      // Compute zeros
      (zeros1,zeros2) := roots(
            poly1,
            poly2,
            n_real);

      // Compute resonance frequencies
      w_min := Modelica.Constants.inf;
      w_max := -Modelica.Constants.inf;
      for i in 1:size(zeros1, 1) loop
        if zeros1[i] <> 0 then
          w := abs(zeros1[i]);
          w_min := min(w_min, w);
          w_max := max(w_max, w);
          w_found := true;
        end if;
      end for;

      for i in 1:size(zeros2, 1) loop
        w :=Modelica.ComplexMath.'abs'(zeros2[i]);
        w_min := min(w_min, w);
        w_max := max(w_max, w);
        w_found := true;
      end for;
    end frequencyRange;

    function frequencyRangeBode
      "Determine min. and max. frequencies for Bode plot"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input ZerosAndPoles tf "ZerosAndPoles transfer function";
      output Modelica.SIunits.AngularVelocity w_min;
      output Modelica.SIunits.AngularVelocity w_max;
    protected
      Real phi_min=Modelica.SIunits.Conversions.from_deg(3);
      Real real_min=1.0e-4;
      Real pi=Modelica.Constants.pi;
      Complex numZeros[:];
      Complex denZeros[:];
      Integer n_num;
      Integer n_den;
      Real w_min1;
      Real w_min2;
      Real w_max1;
      Real w_max2;
    algorithm
      // Compute zeros and poles
      (numZeros,denZeros) := ZerosAndPoles.Analysis.zerosAndPoles(tf);

      // Compute frequencies for numerator
      n_num := size(numZeros, 1);
      if n_num > 0 then
        (w_min1,w_max1) := Modelica_LinearSystems2.Internal.frequencyRangeZeros(
              numZeros,
              phi_min,
              real_min);
      end if;

      // Compute frequencies for denominator
      n_den := size(denZeros, 1);
      if n_den > 0 then
        (w_min2,w_max2) := Modelica_LinearSystems2.Internal.frequencyRangeZeros(
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
    end frequencyRangeBode;

  encapsulated function fromFile_pc
    "Generate a ZerosAndPoles data record by reading the polynomial coefficients from a file (default file name is pc.mat)"
    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input String fileName="pc.mat" "Name of the zeros and poles data file"
      annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));

    protected
    Integer n1n2d1d2[4]=
        ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(fileName);
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

    protected
    Integer n1_2=if n1 > 0 then 1 else 0 "second dimension of n1-matrix";
    Integer n2_2=if n2 > 0 then 2 else 0 "second dimension of n2-matrix";
    Integer d1_2=if d1 > 0 then 1 else 0 "second dimension of d1-matrix";
    Integer d2_2=if d2 > 0 then 2 else 0 "second dimension of d2-matrix";

    Real k = scalar(Streams.readRealMatrix(fileName, "k", 1, 1));
    Real n1Vector[n1] = vector(Streams.readRealMatrix(fileName, "n1", n1, n1_2))
      "Coefficients of first order numenator polynomials";
    Real n2Matrix[n2,n2_2] = Streams.readRealMatrix(fileName, "n2", n2, n2_2)
      "Coefficients of second order denominator polynomials";
    Real d1Vector[d1] = vector(Streams.readRealMatrix(fileName, "d1", d2, d1_2))
      "Coefficients of first order denominator polynomials";
    Real d2Matrix[d2,d2_2] = Streams.readRealMatrix(fileName, "d2", d2, d2_2)
      "Coefficients of second order numenator polynomials";

  algorithm
    zp.k := k;
    zp.n1 := if n1 > 0 then n1Vector else fill(0, 0);
    zp.n2 := if n2 > 0 then n2Matrix else fill(0, 0, 2);
    zp.d1 := if d1 > 0 then d1Vector else fill(0, 0);
    zp.d2 := if d2 > 0 then d2Matrix else fill(0, 0, 2);

  end fromFile_pc;

  encapsulated function fromFile_zp
    "Generate a ZerosAndPoles data record by reading poles and zeros from a file (default file name is zp.mat)"

    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Complex;

    input String fileName = "zp.mat" "Name of the zeros and poles data file"
      annotation (
        Dialog(
          loadSelector(
            filter="MAT files (*.mat);; All files (*.*)",
            caption="state space system data file")));
    protected
    Integer n1n2d1d2[4]=
        ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(fileName);
    Integer n1=n1n2d1d2[1];
    Integer n2=n1n2d1d2[2];
    Integer d1=n1n2d1d2[3];
    Integer d2=n1n2d1d2[4];
    Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
    Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
    output ZerosAndPoles zp(
      n1=fill(0, n1),
      n2=fill(0, n2, 2),
      d1=fill(0, d1),
      d2=fill(0, d2, 2));

    protected
    Integer z_2 = if zSize > 0 then 2 else 0 "second dimension of zeros-matrix";
    Integer p_2 = if pSize > 0 then 2 else 0 "second dimension of poles-matrix";

    Real k = scalar(Streams.readRealMatrix(fileName, "k", 1, 1));
    Real zerosMatrix[zSize,z_2] = Streams.readRealMatrix(fileName, "z", zSize, z_2)
      "Zeros in rows of real parts and imaginary parts";
    Real polesMatrix[pSize,p_2] = Streams.readRealMatrix(fileName, "p", pSize, p_2)
      "Poles in rows of real parts and imaginary parts";
    Complex zeros[:]=if zSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        zerosMatrix[:, 1], zerosMatrix[:, z_2]) else fill(Complex(0), 0);
    Complex poles[:]=if pSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        polesMatrix[:, 1], polesMatrix[:, p_2]) else fill(Complex(0), 0);

  algorithm
    zp := ZerosAndPoles(
        k=k,
        z=zeros,
        p=poles);
  end fromFile_zp;

    function fromRealAndImag
      "Generate a complex vector from a real part vector and imaginary part vector "

      import Complex;

      input Real real[:];
      input Real imag[size(real, 1)];
      output Complex result[size(real, 1)] "Number of real zeros";
    algorithm
      for i in 1:size(real, 1) loop
        result[i].re := real[i];
        result[i].im := imag[i];
      end for;
    end fromRealAndImag;

    encapsulated function isControllableAndObservableSISO
      "To check whether a SISO system is controllable and observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

    input ZerosAndPoles zp;

      output Boolean controllableAndObservable;
    protected
      StateSpace ss=StateSpace(zp);

    algorithm
      controllableAndObservable := StateSpace.Internal.isControllableAndObservableSISO(ss=ss);

    end isControllableAndObservableSISO;

    function isRoot
      "Check if frequency is an element of the complex vector zeros"

      import Complex;
      //import Modelica_LinearSystems2;

      input Complex zeros[:];
      input Complex p;
    //Never used
      input Real eps(min=0) = 0;
      output Boolean result;

    protected
      Integer sz=size(zeros, 1);
      Integer i;

    algorithm
      i := 1;
      result := false;
      while i <= sz and not result loop

    //The implementation of Complex.'==' does not take a third parameter.
        result := p == zeros[i];
        i := i + 1;
      end while;

    end isRoot;

    encapsulated function normalizationFactor
      "Compute correction factor of low pass filter such that amplitude at cut-off frequency is -3db (=10^(-3/20) = 0.70794...)"
      import Modelica;
      import Modelica.Utilities.Streams;

      input Real c1[:]
        "[p] coefficients of denominator polynomials (c1[i}*p + 1)";
      input Real c2[:,2]
        "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
      output Real alpha "Correction factor (replace p by alpha*p)";
    protected
      Real alpha_min;
      Real alpha_max;

      function normalizationResidue "Residue of correction factor computation"
        input Real c1[:]
          "[p] coefficients of denominator polynomials (c1[i]*p + 1)";
        input Real c2[:,2]
          "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
        input Real alpha;
        output Real residue;
      protected
        constant Real beta= 10^(-3/20)
          "Amplitude of -3db required, i.e. -3db = 20*log(beta)";
        Real cc1;
        Real cc2;
        Real p;
        Real alpha2=alpha*alpha;
        Real alpha4=alpha2*alpha2;
        Real A2=1.0;
      algorithm
        assert(size(c1,1) <= 1, "Internal error 2 (should not occur)");
        if size(c1, 1) == 1 then
          cc1 := c1[1]*c1[1];
          p := 1 + cc1*alpha2;
          A2 := A2*p;
        end if;
        for i in 1:size(c2, 1) loop
          cc1 := c2[i, 2]*c2[i, 2] - 2*c2[i, 1];
          cc2 := c2[i, 1]*c2[i, 1];
          p := 1 + cc1*alpha2 + cc2*alpha4;
          A2 := A2*p;
        end for;
        residue := 1/sqrt(A2) - beta;
      end normalizationResidue;

      function findInterval "Find interval for the root"
        input Real c1[:]
          "[p] coefficients of denominator polynomials (a*p + 1)";
        input Real c2[:,2]
          "[p^2, p] coefficients of denominator polynomials (b*p^2 + a*p + 1)";
        output Real alpha_min;
        output Real alpha_max;
      protected
        Real alpha = 1.0;
        Real residue;
      algorithm
        alpha_min :=0;
        residue := normalizationResidue(c1, c2, alpha);
        if residue < 0 then
           alpha_max :=alpha;
        else
           while residue >= 0 loop
              alpha := 1.1*alpha;
              residue := normalizationResidue(c1, c2, alpha);
           end while;
           alpha_max :=alpha;
        end if;
      end findInterval;

    function solveOneNonlinearEquation
        "Solve f(u) = 0; f(u_min) and f(u_max) must have different signs"
        import Modelica.Utilities.Streams.error;

      input Real c1[:]
          "[p] coefficients of denominator polynomials (c1[i]*p + 1)";
      input Real c2[:,2]
          "[p^2, p] coefficients of denominator polynomials (c2[i,1]*p^2 + c2[i,2]*p + 1)";
      input Real u_min "Lower bound of search intervall";
      input Real u_max "Upper bound of search intervall";
      input Real tolerance=100*Modelica.Constants.eps
          "Relative tolerance of solution u";
      output Real u "Value of independent variable so that f(u) = 0";

      protected
      constant Real eps=Modelica.Constants.eps "machine epsilon";
      Real a=u_min "Current best minimum interval value";
      Real b=u_max "Current best maximum interval value";
      Real c "Intermediate point a <= c <= b";
      Real d;
      Real e "b - a";
      Real m;
      Real s;
      Real p;
      Real q;
      Real r;
      Real tol;
      Real fa "= f(a)";
      Real fb "= f(b)";
      Real fc;
      Boolean found=false;
    algorithm
      // Check that f(u_min) and f(u_max) have different sign
      fa := normalizationResidue(c1,c2,u_min);
      fb := normalizationResidue(c1,c2,u_max);
      fc := fb;
      if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
        error(
          "The arguments u_min and u_max to solveOneNonlinearEquation(..)\n" +
          "do not bracket the root of the single non-linear equation:\n" +
          "  u_min  = " + String(u_min) + "\n" + "  u_max  = " + String(u_max)
           + "\n" + "  fa = f(u_min) = " + String(fa) + "\n" +
          "  fb = f(u_max) = " + String(fb) + "\n" +
          "fa and fb must have opposite sign which is not the case");
      end if;

      // Initialize variables
      c := a;
      fc := fa;
      e := b - a;
      d := e;

      // Search loop
      while not found loop
        if abs(fc) < abs(fb) then
          a := b;
          b := c;
          c := a;
          fa := fb;
          fb := fc;
          fc := fa;
        end if;

        tol := 2*eps*abs(b) + tolerance;
        m := (c - b)/2;

        if abs(m) <= tol or fb == 0.0 then
          // root found (interval is small enough)
          found := true;
          u := b;
        else
          // Determine if a bisection is needed
          if abs(e) < tol or abs(fa) <= abs(fb) then
            e := m;
            d := e;
          else
            s := fb/fa;
            if a == c then
              // linear interpolation
              p := 2*m*s;
              q := 1 - s;
            else
              // inverse quadratic interpolation
              q := fa/fc;
              r := fb/fc;
              p := s*(2*m*q*(q - r) - (b - a)*(r - 1));
              q := (q - 1)*(r - 1)*(s - 1);
            end if;

            if p > 0 then
              q := -q;
            else
              p := -p;
            end if;

            s := e;
            e := d;
            if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q) then
              // interpolation successful
              d := p/q;
            else
              // use bi-section
              e := m;
              d := e;
            end if;
          end if;

          // Best guess value is defined as "a"
          a := b;
          fa := fb;
          b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
          fb := normalizationResidue(c1,c2,b);

          if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
            // initialize variables
            c := a;
            fc := fa;
            e := b - a;
            d := e;
          end if;
        end if;
      end while;

      annotation (Documentation(info="<html>

<p>
This function determines the solution of <b>one non-linear algebraic equation</b> &quot;y=f(u)&quot;
in <b>one unknown</b> &quot;u&quot; in a reliable way. It is one of the best numerical
algorithms for this purpose. As input, the nonlinear function f(u)
has to be given, as well as an interval u_min, u_max that
contains the solution, i.e., &quot;f(u_min)&quot; and &quot;f(u_max)&quot; must
have a different sign. If possible, a smaller interval is computed by
inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last 3 points and computing the zero). If this fails,
bisection is used, which always reduces the interval by a factor of 2.
The inverse quadratic interpolation method has superlinear convergence.
This is roughly the same convergence rate as a globally convergent Newton
method, but without the need to compute derivatives of the non-linear
function. The solver function is a direct mapping of the Algol 60 procedure
&quot;zero&quot; to Modelica, from:
</p>

<dl>
<dt> Brent R.P.:</dt>
<dd> <b>Algorithms for Minimization without derivatives</b>.
     Prentice Hall, 1973, pp. 58-59.</dd>
</dl>
</html>"));
    end solveOneNonlinearEquation;

    algorithm
       // Find interval for alpha
       (alpha_min, alpha_max) :=findInterval(c1, c2);

       // Compute alpha, so that abs(G(p)) = -3db
       alpha :=solveOneNonlinearEquation(
        c1,
        c2,
        alpha_min,
        alpha_max);
    end normalizationFactor;

    function numberOfRealPoles "Calculate number of real poles"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input TransferFunction tf "TransferFunction";
    output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(tf.d)));
    algorithm
    end numberOfRealPoles;

    encapsulated function numberOfRealZeros "Calculate number of real zeros"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input Real poly1[:];
      input Real poly2[:,2];
      output Integer result "Number of real zeros";
    protected
      Real D;
    algorithm
      result := size(poly1, 1);
      for i in 1:size(poly2, 1) loop
        D := (poly2[i, 1]/2)^2 - poly2[i, 2];
        if D >= 0 then
            // two real zeros
          result := result + 2;
        end if;
      end for;
    end numberOfRealZeros;

    function numberOfRealZeros2 "Calculate number of real zeros"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input TransferFunction tf "TransferFunction";
      output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(tf.n)));
    algorithm
    end numberOfRealZeros2;

    function numberOfRealZerosAndPoles_pc
      "Generate a zeros and poles data record by reading the polynomial coefficients from a file (default file name is zp.mat)"
      import Modelica;
      import Modelica.Utilities.Streams;

      input String fileName = "pc.mat" "Name of the zeros and poles data file"
        annotation (
          Dialog(
            loadSelector(
              filter="MAT files (*.mat);; All files (*.*)",
              caption="State space system data file")));
      output Integer n1n2d1d2[4];

    protected
      Integer n1Size[2] = Streams.readMatrixSize(fileName, "n1");
      Integer n2Size[2] = Streams.readMatrixSize(fileName, "n2");
      Integer d1Size[2] = Streams.readMatrixSize(fileName, "d1");
      Integer d2Size[2] = Streams.readMatrixSize(fileName, "d2");

    algorithm
      n1n2d1d2[1] := n1Size[1];
      n1n2d1d2[2] := n2Size[1];
      n1n2d1d2[3] := d1Size[1];
      n1n2d1d2[4] := d2Size[1];

    end numberOfRealZerosAndPoles_pc;

    function numberOfRealZerosAndPoles_zp
      "Get the number of first oder polynomials (n1, d1) and second order polynomials (n2, d2) of zeros and poles from zeros and poles written in a MAT-file"

      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.DataDir;
      import Modelica_LinearSystems2.Internal;

      input String fileName=DataDir + "/zp.mat"
        "Name of the zeros and poles data file"      annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
      input Real eps=Modelica.Constants.eps;

      output Integer n1n2d1d2[4];

    protected
      Integer n1;
      Integer d1;

      Integer zSize[2] = Streams.readMatrixSize(fileName, "z");
      Integer pSize[2] = Streams.readMatrixSize(fileName, "p");

      Real zerosMatrix[zSize[1],zSize[2]] = Streams.readRealMatrix(
              fileName,
              "z",
              zSize[1],
              zSize[2]) "zeros in rows of real parts and imaginary parts";
      Real polesMatrix[pSize[1],pSize[2]] = Streams.readRealMatrix(
              fileName,
              "p",
              pSize[1],
              pSize[2]) "poles in rows of real parts and imaginary parts";

    algorithm
      n1 := zSize[1];
      d1 := pSize[1];
      for i in 1:zSize[1] loop
        if abs(zerosMatrix[i, 2]) >= eps then
          n1 := n1 - 1;
        end if;
      end for;
      for i in 1:pSize[1] loop
        if abs(polesMatrix[i, 2]) >= eps then
          d1 := d1 - 1;
        end if;
      end for;
      n1n2d1d2[1] := n1;
      n1n2d1d2[2] := div((zSize[1] - n1), 2);
      n1n2d1d2[3] := d1;
      n1n2d1d2[4] := div((pSize[1] - d1), 2);

    end numberOfRealZerosAndPoles_zp;

    function 'p+a' "Addition of a complex number and a real value"
      import Modelica;
      import Complex;
      import Modelica.ComplexMath.j;

      input Complex p; // "Complex number";
      input Real a "Value of Real variable";
      output Complex c;
    algorithm
      c := p.re + p.im*j +a;
      annotation(Inline=true);
    end 'p+a';

    function 'p^2+k[1]*p+k[2]'
      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica.ComplexMath.j;

      input Complex p;
      input Real k[2];
      output Complex c;
    algorithm
      c := p.re^2 - p.im^2 + k[1]*p.re + k[2]+p.im*(2*p.re + k[1])*j;
      annotation(Inline=true);
    end 'p^2+k[1]*p+k[2]';

    function roots "Determine zeros of factorized polynomial"
      import Modelica;
      import Complex;
      import Modelica.ComplexMath.j;

      input Real poly1[:] "[p^0] coefficients of first order polynomials";
      input Real poly2[:,2] "[p, p^0] coefficients of second order polynomials";
      input Integer n_real
        "Number of real zeros computed with Internal.numberOfRealZeros";
      output Real realZeros[n_real] "All real zeros of poly1 and poly2";
      output Complex complexZeros[:]=fill(Complex(0, 0), integer((size(poly1, 1)
           + 2*size(poly2, 1) - n_real)/2))
        "All complex zeros of poly1 and poly2; for a complex conjugate pair, only one zero is stored";
    protected
      Integer np1=size(poly1, 1);
      Integer np2=size(poly2, 1);
      Real D;
      Real D2;
      Real b;
      Integer j1;
      Integer j2;

    algorithm
      assert(np1 <= n_real, "Size of poly1 = " + String(np1) + " > n_real " +
        " (= " + String(n_real) + ").");
      for i in 1:np1 loop
        realZeros[i] := -poly1[i];
      end for;

      j1 := np1 + 1;
      j2 := 1;
      for i in 1:np2 loop
        b := poly2[i, 1]/2;
        D := b*b - poly2[i, 2];
        D2 := sqrt(abs(D));
        if D >= 0 then
          realZeros[j1] := -b + D2;
          realZeros[j1 + 1] := -b - D2;
          j1 := j1 + 2;
        else
          complexZeros[j2] := -b+D2*j;
          j2 := j2 + 1;
        end if;
      end for;
    end roots;

    function sameMatrixRows "Determine identical rows of a [:,2] matrix"

      input Real M[:,2] "Matrix";
      input Integer startIndex=1 "Start index";
      output Integer endIndex=startIndex
        "startIndex:endIndex are identical rows in M";
    protected
      Integer n=size(M, 1);
      Integer i=startIndex + 1;
      Real v0[size(M, 2)]=M[startIndex, :];
    algorithm
      while i <= n loop
        if M[i, 1] == v0[1] and M[i, 2] == v0[2] then
          endIndex := i;
          i := i + 1;
        else
          i := n + 1;
        end if;
      end while;
    end sameMatrixRows;

    function sameVectorElements "Determine identical elements of a vector"

      input Real v[:] "Vector";
      input Integer startIndex=1 "Start index";
      output Integer endIndex=startIndex
        "startIndex:endIndex are identical elements in v";
    protected
      Integer nv=size(v, 1);
      Integer i=startIndex + 1;
      Real v0=v[startIndex];
    algorithm
      while i <= nv loop
        if v[i] == v0 then
          endIndex := i;
          i := i + 1;
        else
          i := nv + 1;
        end if;
      end while;
    end sameVectorElements;

  encapsulated function scaleFactor1
      "Return scale factor for first order block"
      import Modelica;
    input Real n "(s+n)/(s+d)";
    input Real d "(s+n)/(s+d)";
    input Real small=100*Modelica.Constants.eps;
    output Real k "= n/d, if d,n are not zero, otherwise special cases";
  algorithm
    k := if abs(d) > small  and abs(n) > small then abs(n)/abs(d) else 1;

  end scaleFactor1;

  function scaleFactor2 "Return scale factor for second order block"
      import Modelica;
    input Real n1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
    input Real n2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
    input Real d1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
    input Real d2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
    input Real small=100*Modelica.Constants.eps;
    output Real k "= n2/d2, if d2,n2 are not zero, otherwise special cases";
  algorithm

    k := if abs(n2) > small and abs(d2) > small then d2/n2 else 1;

  end scaleFactor2;

    function secondOrderToString
      "Transform vector of coefficients of second order polynomials to a string representation"
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math;

      input Real c[:,2]=fill(0.0,0,2)
        "[p,p^0] coefficients osecond order polynomials: polynom(p) = p^2 + c[:,1]*p + c[:,2]";
      input Integer significantDigits=6
        "Number of significant digits that are shown";
      input String name="p" "Independent variable name used for printing";
      input Boolean normalized=false
        "= true, the polynomials in the string are represented as p^2/c[:,2] + c[:,1]/c[:,2]*p + 1, provided c[i,2]<>0";
      output String s="";
      output Real gain=1.0
        "If normalized=true, the product(c[i,2]) (for i, where c[i,2]<>0), otherwise gain=1.0";
    protected
      Integer nc=size(c, 1);
      Real cs[nc,2];
      Real cc[nc,2];
      Integer i=1;
      Integer j;
      Integer nj;
      constant Real smallNumber = 100*Modelica.Constants.small;
      constant Real eps = 10*Modelica.Constants.eps;
    algorithm
      // Change coefficients, if normalized output, and sort them
      for i in 1:nc loop
         if Math.isEqual(c[i,2], 0.0, smallNumber) then
            cs[i,1] := c[i,1];
            cs[i,2] := 0.0;
         elseif normalized then
            cs[i,1] := 1/c[i,2];
            cs[i,2] := c[i,1]/c[i,2];
            gain := gain*c[i,2];
         else
            cs[i,1] := c[i,1];
            cs[i,2] := c[i,2];
         end if;
      end for;
      cc :=Modelica.Math.Matrices.sort(cs);

    // Generate string
    if normalized then
      while i <= nc loop
        if i > 1 then
          s := s + "*";
        end if;
        j := sameMatrixRows(cc, i);
        nj := j - i + 1;
        if cc[i, 1] == 0.0 and cc[i, 2] == 0.0 then
          // case p^2
          s := s + name + "^" + String(2*nj);
        else
          // b*p^2 term
          if Math.isEqual(cc[i, 1], 1.0) then
            s := s + "(" + name + "^2";
          elseif Math.isEqual(cc[i, 1], -1.0) then
            s := s + "(-" + name + "^2";
          else
            s := s + "(" + String(cc[i, 1], significantDigits=significantDigits) +
              "*" + name + "^2";
          end if;

          // a*p term
          if cc[i, 2] == 0.0 then
             s := s + " + 1)";
          else
             if Math.isEqual(cc[i, 2], 1.0, eps) then
                s := s + " + ";
             elseif Math.isEqual(cc[i, 2], -1.0, eps) then
                s := s + " - ";
             elseif cc[i, 2] >= 0.0 then
                s := s + " + " + String(cc[i, 2], significantDigits=significantDigits) + "*";
             else
                s := s + " - " + String(-cc[i, 2], significantDigits=significantDigits) + "*";
             end if;
             s := s + name + " + 1)";
          end if;
          if nj > 1 then
            s := s + "^" + String(nj);
          end if;
        end if;
        i := j + 1;
      end while;

    else
      while i <= nc loop
        if i > 1 then
          s := s + "*";
        end if;
        j := sameMatrixRows(cc, i);
        nj := j - i + 1;
        if cc[i, 1] == 0 and cc[i, 2] == 0 then
          // case p^2
          s := s + name + "^" + String(2*nj);
        else
          s := s + "(" + name + "^2";

          // b*p term
          if cc[i, 1] == 1.0 then
            s := s + " + " + name;
          elseif cc[i, 1] == -1.0 then
            s := s + " - " + name;
          elseif cc[i, 1] <> 0.0 then
            if cc[i, 1] > 0 then
              s := s + " + ";
            else
              s := s + " - ";
            end if;
            s := s + String(abs(cc[i, 1]), significantDigits=significantDigits) +
              "*" + name;
          end if;

          // a*p^0 term
          if (cc[i, 2]) > 0.0 then
            s := s + " + " + String(cc[i, 2], significantDigits=significantDigits) + ")";
          elseif (cc[i, 2]) < 0.0 then
            s := s + " - " + String(-cc[i, 2], significantDigits=significantDigits) + ")";
          else
            s := s + ")";
          end if;
          if nj > 1 then
            s := s + "^" + String(nj);
          end if;
        end if;
        i := j + 1;
       end while;
    end if;
    end secondOrderToString;

  end Internal;

  annotation (
    defaultComponentName="filter",
    Documentation(info="<html>
<p>
This record defines a transfer function by its zeros, poles and a gain:
</p>
<blockquote><pre>
         product(p - z[i])
y = k * ------------------- * u
         product(p - n[i])
</pre></blockquote>
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
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Design.filter\">ZerosAndPoles.Design.filter</a>, a
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
<blockquote><pre>
         product(p + n1[i]) * product(p^2 + n2[i,1]*p + n2[i,2])
y = k * ---------------------------------------------------------
         product(p + d1[i]) * product(p^2 + d2[i,1]*p + d2[i,2])
</pre></blockquote>
<p>
Note, the degrees of the numerator and denominator
polynomials are given as:
</p>
<blockquote><pre>
degree of numerator   = size(n1,1) + 2*size(n2,1);
degree of denominator = size(d1,1) + 2*size(d2,1);
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
                         (p+1)
zp = 4 * -------------------------------------
          (p - 1)*(p - (2+j*3))*(p - (2-j*3))
</pre></blockquote>
<p>
with j=sqrt(-1), is defined as
</p>
<blockquote><pre>
<b>import</b> Modelica_LinearSystems2.Math.Complex;
<b>import</b> Modelica_LinearSystems2.ZerosAndPoles;

zp = ZerosAndPoles(z = {Complex(-1,0)},
                   p = {Complex(1,0),
                        Complex(2,3),
                        Complex(2,-3)},
                        k=4);
</pre></blockquote>
</html>"),
    Icon(graphics={
        Rectangle(
          lineColor={160,160,164},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid,
          extent={{-100,-100},{100,100}},
          radius=25.0),
        Text(
          lineColor={255,255,255},
          extent={{-90,-50},{90,50}},
          textString="ZP")}));
end ZerosAndPoles;
