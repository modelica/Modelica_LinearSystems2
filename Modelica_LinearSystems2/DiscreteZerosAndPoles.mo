within Modelica_LinearSystems2;
operator record DiscreteZerosAndPoles
  "Discrete zeros and poles description of a single input, single output system (data + operations)"

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

  encapsulated operator 'constructor'
    "Collection of operators to construct a DiscreteZerosAndPoles data record"
    import Modelica_LinearSystems2;
    import Modelica;

    encapsulated function fromReal
      "Generate a DiscreteZerosAndPoles data record from a real value"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts=0 "Sample time";
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

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = 'constructor'.<b>fromReal</b>(r)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a DiscreteZerosAndPoles record zp from a Real value, i.e. a without dynamics:
</p>
<blockquote><pre>
y = r*u
</pre></blockquote>
<p>
Therefore, the record is defined by
</p>
<blockquote><pre>
  dzp.k = r;
  dzp.n1 = fill(0,1);
  dzp.n2 = fill(0,1,2);
  dzp.d1 = fill(0,1);
  dzp.d2 = fill(0,1,2);
</pre></blockquote>
</html>"));
    end fromReal;

  function fromZerosAndPoles
      "Generate a DiscreteZerosAndPoles data record from a continuous ZerosAndPoles represenation"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input ZerosAndPoles zp "continuous zeros and poles transfer function";
    input Modelica.SIunits.Time Ts "Sample time"
         annotation(Dialog(group="Data used to construct discrete from continuous system"));

      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

    output DiscreteZerosAndPoles dzp;
    protected
    StateSpace ss = StateSpace(zp);
    Modelica_LinearSystems2.DiscreteStateSpace dss=
                             Modelica_LinearSystems2.DiscreteStateSpace(ss,Ts,method);

  algorithm
    dzp := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss);
    annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = 'constructor'.<b>fromZerosAndPoles</b>(zp, Ts, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a DiscreteZerosAndPoles record
</p>
<blockquote><pre>
           product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
dzp = k * ---------------------------------------------------------
           product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre></blockquote>
<p>
from a continuous zeros-and-poles transfer function. The functions applies the discretization methods of
<a href=\"modelica://Modelica_LinearSystems2.DiscreteStateSpace.'constructor'.fromStateSpace\">DiscreteStateSpace.'constructor'.fromStateSpace</a>, i.e,
the continuous zeros-and-poles representation is transformed to continuous state space first and then converted to a discrete state space system. Finally the
discrete zeros-and-poles transfer function is derived from DiscreteStateSpace by
<a href=\"modelica://Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles\">DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  zp=Modelica_LinearSystems2.ZerosAndPoles(      n1={2},  d1={1},      d2=[1,1]);

  dzp=Modelica_LinearSystems2.DiscreteZerosAndPoles.'constructor'.fromZerosAndPoles(zp,0.1)

  //                          (q - 0.818182)*(q^2 + 2*q + 1)
  // dzp = 0.00248841 * ---------------------------------------------
  //                     (q - 0.904762)*(q^2 - 1.89549*q + 0.904988)
</pre></blockquote>
</html>"));
  end fromZerosAndPoles;

  encapsulated function fromPolesAndZeros
      "Generate a DiscreteZerosAndPoles data record from a set of zeros and poles"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Internal;
    import Modelica.Utilities.Streams.print;

    input Complex z[:]=fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
    input Complex p[:]=fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
    input Real k=1.0 "Constant multiplied with transfer function";
    input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
    input String uName="" "input name";
    input String yName="" "output name";
    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[Internal.numberOfRealZeros(z)],
      redeclare Real n2[integer((size(z, 1) - Internal.numberOfRealZeros(z))/2),
        2],
      redeclare Real d1[Internal.numberOfRealZeros(p)],
      redeclare Real d2[integer((size(p, 1) - Internal.numberOfRealZeros(p))/2),
        2]) "ZerosAndPoles transfer functions of the zeros, poles and k";

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
<h4>Syntax</h4>
<blockquote><pre>
dzp = 'constructor'.<b>fromPolesAndZeros</b>(z, p, k, Ts, method)
dzp = 'constructor'.<b>fromPolesAndZeros</b>(z, p, k, Ts, method, uName, yName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a DiscreteZerosAndPoles transfer function from denominator
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
             (q + 0.1)
dzp = 4* -----------------
         (q + 1)*(q - 0.1)
</pre></blockquote>
<p>
is defined as
</p>
<blockquote><pre>
  <b>import</b> Modelica_LinearSystems2.Math.Complex;
  <b>import</b> Modelica_LinearSystems2.DiscreteZerosAndPoles;

  dzp = DiscreteZerosAndPoles(
          z = {Complex(-0.1,0)},
          p = {Complex(-1,0),
               Complex(0.1)},
               k=4);
</pre></blockquote>
</html>"));
  end fromPolesAndZeros;

    function fromDiscreteTransferFunction =
      Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteZerosAndPoles
      "Generate a DiscreteZerosAndPoles data record from a discrete transfer function"
      annotation (Documentation(info="<html>
</html>"));

    encapsulated function fromFactorization
      "Generate a DiscreteZerosAndPoles data record from first and second order polynomials"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

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

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = 'constructor'.<b>fromFactorization</b>(n1, n2, d1, d2, k, Ts, method)
dzp = 'constructor'.<b>fromFactorization</b>(n1, n2, d1, d2, k, Ts, method, uName, yName)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a DiscreteZerosAndPoles transfer function from the real first and scond order polynomials of the numerator and the denominator, respectively.
</p>

<h4>Example</h4>
<blockquote><pre>
              (q + 0.1)
dzp = 4 * -----------------
          (q + 1)*(q - 0.1)
</pre></blockquote>
<p>
is defined as
</p>
<blockquote><pre>
  <b>import</b> Modelica_LinearSystems2.Math.Complex;
  <b>import</b> Modelica_LinearSystems2.DiscreteZerosAndPoles;

  dzp = DiscreteZerosAndPoles(n1={0.1},n2=fill(0,0,2), d1={1,-0.1}, d2=fill(0,0,2), k=4);

which is equal to

  dzp = DiscreteZerosAndPoles(n1={0.1},n2=fill(0,0,2), d1=fill(0,0), d2=[0.9, -0.1], k=4);
</pre></blockquote>
</html>"));
    end fromFactorization;

    annotation (Icon(graphics={
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
  "Contains operators for subtraction of discrete zeros and poles descriptions"
  import Modelica;

  function subtract "Subtract two DiscreteZerosAndPoles (dzp1 - dzp2)"
    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

    protected
    Integer size_z1n1=size(dzp1.n1, 1);
    Integer size_z1d1=size(dzp1.d1, 1);
    Integer size_z1n2=size(dzp1.n2, 1);
    Integer size_z1d2=size(dzp1.d2, 1);
    Integer size_z2n1=size(dzp2.n1, 1);
    Integer size_z2d1=size(dzp2.d1, 1);
    Integer size_z2n2=size(dzp2.n2, 1);
    Integer size_z2d2=size(dzp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output DiscreteZerosAndPoles result "= dzp1-dzp2";
  algorithm
    assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for subtraction with \"-.subtract\".");
    result.Ts := dzp1.Ts;
     if dzp1 == dzp2 then
       result := DiscreteZerosAndPoles(0);
     else
      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,dzp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(1, {1}, dzp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,dzp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(1, {1}, dzp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,dzp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(1, {1}, dzp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,dzp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(1, {1}, dzp1.d2[i, :]));
      end for;

      p3 := dzp1.k*p1 - dzp2.k*p2;
      k := 0;
      for i in size(p3.c, 1):-1:1 loop
        if abs(p3.c[i]) > Modelica.Constants.eps then
          k := p3.c[i];
        end if;
      end for;
      numZeros := Polynomial.roots(p3);
      result := DiscreteZerosAndPoles(numZeros, dummy, k, dzp1.Ts, dzp1.method);

      result.d1 := cat(1, dzp1.d1, dzp2.d1);
      result.d2 := cat(1, dzp1.d2, dzp2.d2);

    end if;
  end subtract;

  function negate "Unary minus (multiply transfer function by -1)"
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

     input DiscreteZerosAndPoles dzp;
     output DiscreteZerosAndPoles result(n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2, k=-dzp.k, Ts=dzp.Ts, method=dzp.method) "= -dzp";
  algorithm
  end negate;
    annotation (Documentation(info="<html>
<p>This package contains operators for subtraction of discrete zeros and poles records. </p>
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

encapsulated operator '*'
  "Contains operators for multiplication of discrete zeros and poles records"
  import Modelica;

function 'dzp*dzp'
      "Multiply two DiscreteZerosAndPoles transfer functions (dzp1 * dzp2)"

      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;

  input DiscreteZerosAndPoles dzp1;
  input DiscreteZerosAndPoles dzp2;
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
  input String uName=dzp1.uName "Input name";
  input String yName=dzp2.yName "Output name";

  output DiscreteZerosAndPoles result "= dzp1 * dzp2";
algorithm
  assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for multiplication with \"*\".");
  result.Ts := dzp1.Ts;
  if dzp1 == DiscreteZerosAndPoles(0) or dzp2 == DiscreteZerosAndPoles(0) then
    result := DiscreteZerosAndPoles(0);
  else
    result.n1 := cat(
      1,
      dzp1.n1,
      dzp2.n1);
    result.n2 := cat(
      1,
      dzp1.n2,
      dzp2.n2);
    result.d1 := cat(
      1,
      dzp1.d1,
      dzp2.d1);
    result.d2 := cat(
      1,
      dzp1.d2,
      dzp2.d2);
    result.k := dzp1.k*dzp2.k;
  end if;
    result.Ts := dzp1.Ts;
    result.method := dzp1.method;
    result.uName := uName;
    result.yName := yName;

end 'dzp*dzp';

function 'r*dzp'
      "Multiply a real number with a discrete DiscreteZerosAndPoles transfer function  (r * dzp2)"

      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;

  input Real r;
  input DiscreteZerosAndPoles dzp;
  input String uName=dzp.uName "Input name";
  input String yName=dzp.yName "Output name";

  output DiscreteZerosAndPoles result "= r * dzp1";
algorithm
  if r == 0 or dzp == DiscreteZerosAndPoles(0) then
    result := DiscreteZerosAndPoles(0);
  else
    result := dzp;
    result.k := r*dzp.k;
  end if;
    result.Ts := dzp.Ts;
    result.method := dzp.method;
    result.uName := uName;
    result.yName := yName;

end 'r*dzp';
    annotation (Documentation(info="<html>
<p>This package contains operators for multiplication of discrete zeros and poles records. </p>
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
          radius=25.0),
        Line(
          points={{-36,36},{36,-36}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
            points={{0,50},{0,-50}},
            color={0,0,0},
            smooth=Smooth.None),
        Line(
          points={{36,36},{-36,-36}},
          color={0,0,0},
          smooth=Smooth.None)}));
end '*';

  encapsulated operator function '+'
    "Addition of to discrete transfer functions dzp1 + dzp2, i.e. parallel connection of two transfer functions (= inputs are the same, outputs of the two systems are added)"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

  protected
    Integer size_z1n1=size(dzp1.n1, 1);
    Integer size_z1d1=size(dzp1.d1, 1);
    Integer size_z1n2=size(dzp1.n2, 1);
    Integer size_z1d2=size(dzp1.d2, 1);
    Integer size_z2n1=size(dzp2.n1, 1);
    Integer size_z2d1=size(dzp2.d1, 1);
    Integer size_z2n2=size(dzp2.n2, 1);
    Integer size_z2d2=size(dzp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output DiscreteZerosAndPoles result "= dzp1+dzp2";

  algorithm
    if max({size(dzp1.n1,1),size(dzp1.n2,1),size(dzp1.d1,1),size(dzp1.d2,1)})>0 and max({size(dzp2.n1,1),size(dzp2.n2,1),size(dzp2.d1,1),size(dzp2.d2,1)})>0 then
      assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for addition with \"+\".");
    end if;
    result.Ts := dzp1.Ts;

    if dzp1 == -dzp2 then
      result := DiscreteZerosAndPoles(0,Ts = dzp1.Ts);
    else

      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,dzp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          dzp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,dzp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          dzp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,dzp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          dzp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,dzp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          dzp1.d2[i, :]));
      end for;

      p3 := dzp1.k*p1 + dzp2.k*p2;
      k := p3.c[1];
      numZeros := Polynomial.roots(p3);
      result := DiscreteZerosAndPoles(numZeros, dummy, k, dzp1.Ts, dzp1.method);

      result.d1 := cat(
        1,
        dzp1.d1,
        dzp2.d1);
      result.d2 := cat(
        1,
        dzp1.d2,
        dzp2.d2);
    end if;

  end '+';

  encapsulated operator function '/'
    "Divide two discrete transfer functions in zeros and poles representation (dzp1 / dzp2)"
    import Modelica;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

    output DiscreteZerosAndPoles result "= dzp1/dzp2";

  protected
    Boolean dzp1IsReal=DiscreteZerosAndPoles.Internal.isReal(dzp1);
    Boolean dzp2IsReal=DiscreteZerosAndPoles.Internal.isReal(dzp2);

  algorithm
    assert(abs(dzp2.k) > 100*Modelica.Constants.small, "dzp2 in operator \"Modelica_LinearSystems2.TransferFunction.'/'()\" may not be zero");
  //   if max({size(dzp1.n1, 1),size(dzp1.n2, 1),size(dzp1.d1, 1),size(dzp1.d2, 1)}) >0 and
  //     max({size(dzp2.n1, 1),size(dzp2.n2, 1),size(dzp2.d1, 1),size(dzp2.d2, 1)}) > 0 then
    if not dzp1IsReal and not dzp2IsReal then
        assert(abs(dzp1.Ts - dzp2.Ts) <= Modelica.Constants.eps, "Two discrete zeros-and-poles systems must have the same sample time Ts for division with \"/\".");
    end if;
    result.Ts := if dzp1IsReal then dzp2.Ts else dzp1.Ts;
    result.method := if dzp1IsReal then dzp2.method else dzp1.method;

    if dzp1 == DiscreteZerosAndPoles(0) then
      result := DiscreteZerosAndPoles(0);
    else
      result.n1 := cat(
        1,
        dzp1.n1,
        dzp2.d1);
      result.n2 := cat(
        1,
        dzp1.n2,
        dzp2.d2);
      result.d1 := cat(
        1,
        dzp1.d1,
        dzp2.n1);
      result.d2 := cat(
        1,
        dzp1.d2,
        dzp2.n2);
      result.k := dzp1.k/dzp2.k;
    end if;
  end '/';

  encapsulated operator function '^'
    "Integer power of DiscreteZerosAndPoles (dzp^k)"

    import Modelica;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp;
    input Integer k;

    output DiscreteZerosAndPoles result(
      redeclare Real n1[k*size(dzp.n1, 1)],
      redeclare Real n2[k*size(dzp.n2, 1),2],
      redeclare Real d1[k*size(dzp.d1, 1)],
      redeclare Real d2[k*size(dzp.d2, 1),2]) "= dzp^k";
  protected
    Integer size_n1=size(dzp.n1, 1);
    Integer size_d1=size(dzp.d1, 1);
    Integer size_n2=size(dzp.n2, 1);
    Integer size_d2=size(dzp.d2, 1);
  algorithm
    result.Ts := dzp.Ts;
    result.method := dzp.method;
    for i in 1:k loop
      result.n1[(i - 1)*size_n1 + 1:i*size_n1] := dzp.n1;
      result.d1[(i - 1)*size_d1 + 1:i*size_d1] := dzp.d1;
      result.n2[(i - 1)*size_n2 + 1:i*size_n2, :] := dzp.n2;
      result.d2[(i - 1)*size_d2 + 1:i*size_d2, :] := dzp.d2;
    end for;
    result.k := dzp.k^k;

  end '^';

encapsulated operator function '=='
    "Check whether two discreteZerosAndPoles transfer functions are identical"
    import Modelica;
    import Modelica.Math;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;
    input Real eps(min=0) = 0
      "Two numbers n1 and n2 are identical if abs(n1-n2) <= eps";

    output Boolean result "= dzp1 == dzp2";
algorithm
    result := Math.Vectors.isEqual(dzp1.n1,dzp2.n1,eps) and Math.Vectors.isEqual(dzp1.d1,dzp2.d1,eps) and Math.Matrices.isEqual(dzp1.n2,dzp2.n2,eps) and Math.Matrices.isEqual(dzp1.d2,dzp2.d2,eps) and (dzp1.k==dzp2.k) and abs(dzp1.Ts-dzp2.Ts)<=eps;

end '==';

  encapsulated operator function 'String'
    "Transform DiscreteZerosAndPoles transfer function into a String representation"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Types.Method;

      input DiscreteZerosAndPoles dzp
      "DiscreteZerosAndPoles transfer function to be transformed in a String representation";
      input Integer significantDigits=6
      "Number of significant digits that are shown";
      input String name="q" "Independent variable name used for printing";
      output String s="";
  protected
      Integer num_order=size(dzp.n1, 1) + 2*size(dzp.n2, 1);
      Integer den_order=size(dzp.d1, 1) + 2*size(dzp.d2, 1);
  algorithm
      if num_order == 0 and den_order == 0 then
        s := String(dzp.k);

      else
         // construct string for multiplicative factor
        if dzp.k <> 1.0 or dzp.k == 1.0 and num_order == 0 then
          s := String(dzp.k);
          if num_order <> 0 then
            s := s + "*";
          end if;
        end if;

        if num_order <> 0 then
            // construct numerator string
          s := s + ZerosAndPoles.Internal.firstOrderToString(
                dzp.n1,
                significantDigits,
                name);
          if size(dzp.n2, 1) <> 0 then
            s := if size(dzp.n1, 1) > 0 then s + "*" +
              ZerosAndPoles.Internal.secondOrderToString(
                  dzp.n2,
                  significantDigits,
                  name) else s + ZerosAndPoles.Internal.secondOrderToString(
                  dzp.n2,
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
          s := s + ZerosAndPoles.Internal.firstOrderToString(
                dzp.d1,
                significantDigits,
                name);
          if size(dzp.d2, 1) <> 0 then
            if size(dzp.d1, 1) > 0 then
              s := s + "*";
            end if;
            s := s + ZerosAndPoles.Internal.secondOrderToString(
                  dzp.d2,
                  significantDigits,
                  name);
          end if;
          if den_order > 1 then
            s := s + " )";
          end if;
        end if;
          end if;
          s := s +"\n    Ts = " + String(dzp.Ts) + "\n    method ="+ Modelica_LinearSystems2.Internal.methodString(dzp.method);
    //    end toString;
  end 'String';

  encapsulated function q
    "Generate the DiscreteZerosAndPoles transfer function q"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input Modelica.SIunits.Time Ts=0;
    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[1],
      redeclare Real n2[0,2],
      redeclare Real d1[0],
      redeclare Real d2[0,2]);
  algorithm
    dzp.n1[1] := 0;
    dzp.Ts := Ts;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
q = DiscreteZerosAndPoles.<b>q</b>()
</pre></blockquote>

<h4>Description</h4>
<p>
Generate the complex Laplace variable q=exp(s*T) as a discrete zeros and poles transfer function. It can be used for generating like
</p>
<blockquote><pre>
DiscreteZerosAndPoles dzp = q/(q^2 + q + 1)/(q + 1)
</pre></blockquote>
</html>"));
  end q;

  encapsulated package Analysis
    "Package of functions to analyse discrete zeros-and-poles description represented by a DiscreteZerosAndPoles record"
    import Modelica;
    extends Modelica.Icons.Package;

  encapsulated function timeResponse
      "Calculate the time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;  // Input/Output declarations of discrete time response functions
      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step;
      input Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

    protected
      DiscreteStateSpace dss=DiscreteStateSpace(dzp);
      Real tSpanVar;

  algorithm
    // set sample time
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) := DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=response,
          x0=x0);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x)=DiscreteZerosAndPoles.Analysis.<b>timeResponse</b>(dzp, tSpan, responseType, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the DiscreteZerosAndPoles record is transformed into discrete
state space representation which is given to DiscreteStateSpace.Analysis.timeResponse
to calculate the time response of the state space system. The type of the time
response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step    &quot;Step response&quot;,
Ramp    &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;
</pre></blockquote>
<p>
Starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and x are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  q=Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=1/(q^2 + q + 1)
  dzp.Ts=0.1;

  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[2]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.timeResponse(dzp,tSpan, response,x0);
//  y[:,1,1]={0, 0, 1, 0, 0}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0, 3.0, 0, 0}
</pre></blockquote>
</html>"));
  end timeResponse;

  encapsulated function impulseResponse
      "Calculate the impulse time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
          x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

      annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x)=DiscreteZerosAndPoles.Analysis.<b>impulseResponse</b>(dzp, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>impulseResponse</b> calculates the time response of
a DiscreteZerosAndPoles transfer function with impulse imput.
After transforming the DiscreteZerosAndPoles representation to
DiscreteStateSpace the values of the impulse response are calculated
starting from
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0 for each time step t=k*dt.
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.impulseResponse(dzp, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.timeResponse(dzp, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=dzp=1/(q^2 + q + 1)
  dzp.Ts=0.1;

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.impulseResponse(dzp,tSpan);
//  y[:,1,1]={0, 0, 1, -1.0, 2.96059473233375E-016}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0, 3.0, -3.0, 8.88178419700125E-016}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.timeResponse\">DiscreteZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end impulseResponse;

  encapsulated function stepResponse
      "Calculate the step time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
    Real tSpanVar;
  algorithm

  // set simulation time span
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(DiscreteStateSpace(dzp));
    else
      tSpanVar := tSpan;
    end if;

    (y,t,x_discrete) :=DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
          x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x)=DiscreteZerosAndPoles.Analysis.<b>stepResponse</b>(dzp, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>stepResponse</b> calculates the step response of a DiscreteZerosAndPoles transfer function.
The system is transformed to a DiscreteStateSapce system and starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.stepResponse(dzp, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.timeResponse(dzp, tSpan, response=Types.TimeResponse.Step, x0=fill(0,DiscreteZerosAndPoles.Analysis.denominatorDegree(zp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=dzp=1/(q^2 + q + 1)
  dzp.Ts=0.1;

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.stepResponse(dzp,tSpan);
//  y[:,1,1] = {0, 0, 1, 0, 2.96059473233375E-016}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0, 0, 3.0, 0, 8.88178419700125E-016}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.timeResponse\">DiscreteZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end stepResponse;

  encapsulated function rampResponse
      "Calculate the ramp time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
          x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

      annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x)=DiscreteZerosAndPoles.Analysis.<b>rampResponse</b>(ss, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rampResponse</b> calculates the time response of a DiscreteZerosAndPoles transfer function for ramp imput u = t.
The system is transformed to DiscreteStateSpace state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.rampResponse(dzp, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.timeResponse(dzp, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=dzp=1/(q^2 + q + 1)
  dzp.Ts=0.1;

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.rampResponse(dzp,tSpan);
//  y[:,1,1] = {0, 0, 0, 0.1, 0.1}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1 = {0, 0, 0, 0.3, 0.3}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.timeResponse\">DiscreteZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end rampResponse;

  encapsulated function initialResponse
      "Calculate the initial time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input Real x0[:]=fill(0, 0) "Initial state vector";

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
          x0=x0);

      annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x)=DiscreteZerosAndPoles.Analysis.<b>initialResponse</b>(zp, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>initialResponse</b> calculates the time response of a state space system for given initial condition and zero inputs.
The system is transformed a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.initialResponse(x0, dzp, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteZerosAndPoles.Analysis.timeResponse(dzp, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=dzp=1/(p^2 + p + 1)
  dzp.Ts=0.1;
  Real tSpan= 0.4;
  Real x0[2] = {1,1};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.initialResponse(x0,dzp,tSpan);
//  y[:,1,1 = {0.333333333333, 0.333333333333, -0.666666666667, 0.333333333333,  0.333333333333}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {1, 1, -2.0, 1.0, 1}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.timeResponse\">DiscreteZerosAndPoles.Analysis.timeResponse</a>
</p>
</html>"));
  end initialResponse;

  encapsulated function denominatorDegree
      "Return denominator degree of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Integer result;
  algorithm
    result := size(dzp.d1, 1) + 2*size(dzp.d2, 1);
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = DiscreteZerosAndPoles.Analysis.<b>denominatorDegree</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of
the denominator polynomial constituted by the first and second order
polynomials of the DiscreteZeroAndPoles denominator.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles p = Modelica_LinearSystems2.DiscreteZerosAndPoles.p();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q+0.1)/(p^2 + 0.8*p+ 0.2);

  Real dDegree;

<b>algorithm</b>
  dDegree := DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp);
//  dDegree = 2
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.numeratorDegree\">DiscreteZerosAndPoles.Analysis.numeratorDegree</a>
</p>
</html>"));
  end denominatorDegree;

  encapsulated function numeratorDegree
      "Return numerator degree of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Integer result;
  algorithm
    result := size(dzp.n1, 1) + 2*size(dzp.n2, 1);
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = DiscreteZerosAndPoles.Analysis.<b>numeratorDegree</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>numeratorDegree</b> calculates the degree
of the numerator polynomial constituted by the first and second
order polynomials of the DiscreteZeroAndPoles numerator.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q+0.1)/(p^2 + 0.8*p - 0.2);

  Real nDegree;

<b>algorithm</b>
  nDegree := DiscreteZerosAndPoles.Analysis.numeratorDegree(dzp);
//  nDegree = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.denominatorDegree\">DiscreteZerosAndPoles.Analysis.denominatorDegree</a>.
</p>
</html>"));
  end numeratorDegree;

  encapsulated function evaluate
    "Evaluate a DiscreteZerosAndPoles transfer function at a given value of q"
    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.ZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    input Complex q=Complex(0) "Complex value q where dzp is evaluated";
    input Real den_min=0 "|denominator(p)| is limited by den_min";
    output Complex y "= zp(p)";
    protected
    Complex j = Modelica.ComplexMath.j;
    Complex num;
    Complex den;
    Real abs_den;
  algorithm
    // Build numerator
    num := dzp.k+0*j;
    for i in 1:size(dzp.n1, 1) loop
       num := num*ZerosAndPoles.Internal.'p+a'(q, dzp.n1[i]);
    end for;
    for i in 1:size(dzp.n2, 1) loop
       num := num*ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(q, dzp.n2[i, :]);
    end for;

    // Build denominator
    den := 1+0*j;
    for i in 1:size(dzp.d1, 1) loop
      den := den*ZerosAndPoles.Internal.'p+a'(q, dzp.d1[i]);
    end for;
    for i in 1:size(dzp.d2, 1) loop
      den := den*ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(q, dzp.d2[i, :]);
    end for;

    // Build value of transfer function
    abs_den :=Modelica.ComplexMath.'abs'(den);
    den := if abs_den >= den_min then den else -abs_den+0*j;
    y := num/den;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = DiscreteZerosAndPoles.Analysis.<b>evaluate</b>(dzp,q)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the DiscreteZerosAndPoles
transfer function at a given (complex) value of q.
The transfer function G(z)=N(q)/D(q) is evaluated by calculating the
numerator polynomial N(z) and the denominator polynomial D(q).
</p>

<h4>Example</h4>
<blockquote><pre>
  Complex j = Modelica_LinearSystems2.Math.Complex.j();
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q+1)/(q^2+q+1);

  Complex result;

<b>algorithm</b>
  result := Modelica_LinearSystems2.DiscreteZerosAndPoles.Analysis.evaluate(dzp, j+1);
//  result = 0.538462 - 0.307692j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>
</html>"));
  end evaluate;

  end Analysis;

  encapsulated package Design
    "Package of functions to design discrete zeros-and-poles controllers and observers"
    import Modelica;
    extends Modelica.Icons.Package;

  end Design;

  encapsulated package Plot
    "Package of functions to plot discrete zeros and poles description responses"
    import Modelica;
    extends Modelica.Icons.Package;

  encapsulated function bode
    "Plot discrete zeros a-and-poles transfer function as bode plot"
    import Modelica;
    import Modelica.Utilities.Strings;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.Internal;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Plot;
    import SI = Modelica.SIunits;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles function to be plotted";
    input Integer nPoints(min=2) = 200 "Number of points";
    input Boolean autoRange=true
        "True, if abszissa range is automatically determined";
    input SI.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false"
                                                      annotation(Dialog(enable=not autoRange));
    input SI.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                                annotation(Dialog(enable=not autoRange));

    input Boolean magnitude=true "= true, to plot the magnitude of dzp"
                                                                       annotation(choices(checkBox=true));
    input Boolean phase=true "= true, to plot the pase of dzp" annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot: "
           + String(dzp)));

    input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)"
                                                                           annotation(choices(checkBox=true));
    input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))"
                                                                              annotation(choices(checkBox=true),Dialog(enable=magnitude));

    protected
    SI.AngularVelocity w[nPoints];
    Complex z[nPoints];
    SI.Frequency f[nPoints];
    SI.Conversions.NonSIunits.Angle_deg phi[nPoints];
    Real A[nPoints];
    Boolean OK;
    Complex c;
    SI.Angle phi_old;
    Complex numZeros[:];
    Complex denZeros[:];
    Complex numZerosZ[:];
    Complex denZerosZ[:];
    ZerosAndPoles zp=ZerosAndPoles(k=dzp.k, n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2);

    Plot.Records.Curve curves[2];
    Integer i;
    Plot.Records.Diagram diagram2[2];

  algorithm
          // Determine frequency vector f
    if autoRange then
       (numZerosZ,denZerosZ) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
    else
      numZerosZ := fill(Complex(0), 0);
      denZerosZ := fill(Complex(0), 0);
    end if;

    numZeros := fill(Complex(0),0);
    denZeros := fill(Complex(0),size(denZerosZ,1));
    for i in 1:size(denZerosZ,1) loop
      denZeros[i] :=Modelica.ComplexMath.log(denZerosZ[i])/dzp.Ts;
    end for;

    f := Internal.frequencyVector(
          nPoints,
          autoRange,
          f_min,
          f_max,
          numZeros,
          denZeros);

    // Compute magnitude/phase at the frequency points
    phi_old := 0.0;
    for i in 1:nPoints loop
      w[i] := SI.Conversions.from_Hz(f[i]);
      z[i] :=Modelica.ComplexMath.exp(Complex(0, w[i]*dzp.Ts));
      c := ZerosAndPoles.Analysis.evaluate(
            zp,
            z[i],
            1e-10);
      A[i] :=Modelica.ComplexMath.'abs'(c);
      phi_old :=Modelica.ComplexMath.arg(c, phi_old);
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

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>bode</b>(dzp)
   or
DiscreteZerosAndPoles.Plot.<b>bode</b>(
  dzp,
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
This function plots the bode-diagram of a DiscreteZerosAndPoles transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/( (q - 0.8)*(q - 0.75) )

<b>algorithm</b>
  Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.plotBode(dzp)
//  gives:
</pre></blockquote>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/dBbodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/dBodePhase.png\">
</p>
</html>"));
  end bode;

  encapsulated function timeResponse
      "Plot the time response of a system represented by a discrete zeros-and-poles transfer function. The response type is selectable"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;
  //  input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "Type of time response";
    input Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Time response of  dzp = "
           + String(dzp)));

    protected
    Plot.Records.Curve curve;
    Plot.Records.Diagram diagram2;
    Real y[:,1,1] "Output response";
    Real t[:] "Time vector: (number of samples)";

    Real yy[:] "Output response";
    Real tt[:] "Time vector: (number of samples)";

  algorithm
    (y,t) := DiscreteZerosAndPoles.Analysis.timeResponse(
      dzp,
      tSpan,
      response,
      x0);

    tt := fill(0,2*size(t,1)-1);
    yy := fill(0,2*size(t,1)-1);

    for i in 1:size(t,1)-1 loop
      tt[2*i-1] := t[i];
      tt[2*i] := t[i+1];
      yy[2*i-1] := y[i,1,1];
      yy[2*i] := y[i,1,1];
    end for;
    tt[size(tt,1)] := t[size(t,1)];
    yy[size(tt,1)] := y[size(t,1),1,1];

    curve := Plot.Records.Curve(
      x=tt,
      y=yy,
      legend="y",
      autoLine=true);
    diagram2 := defaultDiagram;
    diagram2.curve := {curve};

    Plot.diagram(diagram2, device);
    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>timeResponse</b>(dzp);
   or
DiscreteZerosAndPoles.Plot.<b>timeResponse</b>(
  dzp,
  tSpan,
  response,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the time response of a discrete zeros and poles transfer
function. The character of the time response if defined by the input
<a href=\"modelica://Modelica_LinearSystems2.Types.TimeResponse\">response</a>,
i.e. Impulse, Step, Ramp, or Initial.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/( (q - 0.8)*(q - 0.75) )
  dzp.Ts = 0.1;

  Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
   Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(dzp, tSpan=3, response=response)
//  gives:
</pre></blockquote>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/DiscreteZerosAndPoles/timeResponseDZP.png\">
</p>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end timeResponse;

  encapsulated function impulse
      "Impulse response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteZerosAndPoles dzp "zeros-and-poles transfer function";
    input Real tSpan=0 "Simulation time span [s]";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
      defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Impulse response of  zp = "+String(dzp)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "Type of time response";
    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
      DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

  algorithm
  // set sample time
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
      dzp=dzp,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>impulse</b>(dzp)
   or
DiscreteZerosAndPoles.Plot.<b>impulse</b>(
  dzp,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the impulse response of a discrete zeros-and-poles transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/( (q - 0.8)*(q - 0.75) )
  dzp.Ts = 0.1;

<b>algorithm</b>
  Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse(dzp, tSpan=2)
//  gives:
</pre></blockquote>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/DiscreteZerosAndPoles/impulseResponseDZP.png\">
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end impulse;

  encapsulated function step
      "Step response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  dzp = "
           + String(dzp)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";
  algorithm
   DiscreteZerosAndPoles.Plot.timeResponse(
      dzp=dzp,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

  equation

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>step</b>(dzp)
   or
DiscreteZerosAndPoles.Plot.<b>step</b>(
  dzp,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the step response of a transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/((q - 0.8)*(q - 0.75));
  dzp.Ts = 0.1;

<b>algorithm</b>
  Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp, tSpan=3)
  // gives:
</pre></blockquote>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/DiscreteZerosAndPoles/stepResponseDZP.png\">
</p>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end step;

  encapsulated function ramp
      "Ramp response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  dzp = "
           + String(dzp)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";

    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";
  algorithm
   Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
      dzp=dzp,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>ramp</b>(dzp)
   or
DiscreteZerosAndPoles.Plot.<b>ramp</b>(
  dzp,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the ramp response of a zeros-and-poles transfer function. It is based on
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles p = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/( (q - 0.8)*(q - 0.75) )
  dzp.Ts = 0.1;

<b>algorithm</b>
  Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp(dzp)
  //  gives:
</pre></blockquote>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/DiscreteZerosAndPoles/rampResponseDZP.png\">
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end ramp;

  encapsulated function initialResponse
      "Initial condition response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";
    input Real y0 "Initial output (for initial condition plot)";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Initial response of  dzp = "
           + String(dzp) + "  with y0 = " + String(y0)));

    protected
    Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(dzp);
    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=
        Modelica.Math.Matrices.equalityLeastSquares(
        dss.A,
        fill(0, size(dss.B, 1)),
        dss.C,
        vector(y0)) "Initial state vector (for initial condition plot)";
  algorithm
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
          dzp=dzp,
          tSpan=tSpan,
          response=response,
          x0=x0,
          defaultDiagram=defaultDiagram,
          device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteZerosAndPoles.Plot.<b>initialResponse</b>(zp)
   or
DiscreteZerosAndPoles.Plot.<b>initialResponse</b>(
  zp,
  tSpan,
  y0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the initial response, i.e. the zeros input response of a zeros and poles transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q^2 - 1.5*q + 0.6)/( (q - 0.8)*(q - 0.75) )
  Real y0=1;
  dzp.Ts = 0.02;

  Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.initialResponseZP(dzp, y0=y0, tSpan=1)
  // gives:
</pre></blockquote>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/DiscreteZerosAndPoles/initialResponseDZP.png\">
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp\">ramp</a>
</p>
</html>"));
  end initialResponse;

  end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of DiscreteZerosAndPoles data record"
    import Modelica_LinearSystems2;
    import Modelica;
    extends Modelica.Icons.Package;
  function toDiscreteTransferFunction
    "Generate a DiscreteTransferFunction object from a DiscreteZerosAndPoles object"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.ZerosAndPoles;

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
<blockquote><pre>
dtf = DiscreteZerosAndPoles.Conversion.<b>toDiscreteTransferFunction</b>(dzp)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a DiscreteTransferFunction record
</p>
<blockquote><pre>
       n(z)     b0 + b1*z + ... + bn*z^n
dtf = ------ = --------------------------
       d(z)     a0 + a1*z + ... + an*z^n
</pre></blockquote>
<p>
from a DiscreteZerosAndPoles record representated by first and second order numerator and denominator polynomials. The poles and zeros and the gain <tt>k</tt> are computed (<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs in the DiscreteTransferFunction constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp = 1/(q + 3)/(q + 1)

<b>algorithm</b>
  dtf:=Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp);
//  dtf = 1/(z^2 + 4*z + 3)
</pre></blockquote>
</html>"));
  end toDiscreteTransferFunction;

  function toDiscreteStateSpace
    "Transform a DiscreteZerosAndPoles object into a DiscreteStateSpace object"
   //encapsulated function fromZerosAndPoles
    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "ZerosAndPoles transfer function of a system";
    output DiscreteStateSpace dss(
      redeclare Real A[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp),
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)],
      redeclare Real B[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp),1],
      redeclare Real C[1,DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)],
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
    Real aa "a2 + a1 +1";
    Real bb "b2 + b1 +1";
    Integer nx=max(DiscreteZerosAndPoles.Analysis.numeratorDegree(dzp),DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp));
    Integer n_num1=size(dzp.n1, 1);
    Integer n_num2=size(dzp.n2, 1);
    Integer n_den1=size(dzp.d1, 1);
    Integer n_den2=size(dzp.d2, 1);
    Integer n_num=n_num1 + 2*n_num2;
    Integer n_den=n_den1 + 2*n_den2;

    Integer i_d=if n_num2 > n_den2 then 2*(n_num2 - n_den2) + 1 else 1;
    Integer i_k=if n_num2 > n_den2 then n_den2 - (n_num2 - n_den2) else n_den2;
    Integer i;
    Integer ili=if max(n_den2, n_num2) > 0 then i_d else  max(2, i_d);
    Real num[nx,2]=[dzp.n2; [dzp.n1,zeros(n_num1)]; zeros(max(0,nx - n_num2 - n_num1), 2)]
        "Numerator matrix, in order that indices are defined in all situations in all if clauses";
    Real den[nx,2]=[dzp.d2; [dzp.d1,zeros(n_den1)]; zeros(max(0,nx - n_den2 - n_den1), 2)]
        "Denominator matrix, in order that indices are defined in all situations in all if clauses";
    Real k[i_k + n_den1](each fixed=false)
        "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
    Real k_total;
    Real eps=1e-5;//100*Modelica.Constants.eps;

    Boolean dZero=true;

  algorithm
    assert(n_num <= n_den,
      "DiscreteZerosAndPoles transfer function is not proper as required from StateSpace system:\n"
       + "  numerator degree (= " + String(n_num) +
      ") <= denominator degree (= " + String(n_den) + ") required.");

    if n_den > 0 then
      for i in 1:max(n_den2, n_num2) loop
        // State space systems of order 2
        if i <= n_den2 then
          if i <= n_num2 then
              // State space system in form (1)
            k[i] := DiscreteZerosAndPoles.Internal.scaleFactor2(
                num[i, 1],
                num[i, 2],
                den[i, 1],
                den[i, 2],eps);
          elseif 2*(i - n_num2) <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
            k[i] := DiscreteZerosAndPoles.Internal.scaleFactor2(
                num[2*(i - n_num2)-1, 1] + num[2*(i - n_num2), 1],
                num[2*(i - n_num2)-1, 1]*num[2*(i - n_num2), 1],
                den[i, 1],
                den[i, 2],eps);
          elseif  2*(i-n_num2) -1== n_num1 then
              // State space system in form (2) with 1 first order numerator polynomial
            k[i] := DiscreteZerosAndPoles.Internal.scaleFactor2(
                0,
                num[2*i-n_num2-1, 1],
                den[i, 1],
                den[i, 2],eps);
           else
              // State space system in form (3)
            k[i] := DiscreteZerosAndPoles.Internal.scaleFactor2(
                0,
                0,
                den[i, 1],
                den[i, 2],eps);
          end if;
        else
           // State space system in form (1) with 2 first order denominator polynomials

          k[i] := DiscreteZerosAndPoles.Internal.scaleFactor2(
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
          k[i_k + i] := DiscreteZerosAndPoles.Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
            n_num2) + i), 1], den[n_den2 + i, 1],eps);
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
           // State space system in form (4)
          k[i_k + i] := DiscreteZerosAndPoles.Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
            1], den[n_den2 + i, 1],eps);
        else
           // State space system in form (5)
          k[i_k + i] := DiscreteZerosAndPoles.Internal.scaleFactor1(0, den[n_den2 + i, 1],eps);
        end if;
      end for;

      k_total := dzp.k/product(k);

      dss.A := zeros(nx, nx);
      dss.B := zeros(nx, 1);
      dss.C := zeros(1, nx);
      dss.D := zeros(1, 1);

   // Calculation of matrices A, B, C, D
   //first elements of A, B, C and D

      if max(n_den2, n_num2) > 0 then

        A[1, :] := {0,1};
        B[1, 1] := 0;
            // Construct state space systems of order 2
        if 1 <= n_den2 then
          aa := num[1, 2] + num[1, 1] +1;
          bb := den[1, 2] + den[1, 1] +1;
          A[2, :] := {-den[1, 2],-den[1, 1]};
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;
          if 1 <= n_num2 then
                 // State space system in form (1)
            C := if abs(bb)>eps and abs(aa)>eps then k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/bb else k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;

          elseif 1 - n_num2 + 1 <= n_num1 then
                 // State space system in form (1) with 2 first order numerator polynomials
            aa := num[1, 1]*num[2, 1] + num[1, 1] + num[2, 1];
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/bb else k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;
          elseif 1 - n_num2 == n_num1 then
                 // State space system in form (2) with 1 first order numerator polynomial
            aa := num[1, 1] + 1;
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps  and abs(aa)>eps then k[1]*[num[1, 1],1]/bb else k[1]*[num[1, 1],1];
            D := [0];
            dZero := dZero and true;
          else
                 // State space system in form (3)
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps  and abs(aa)>eps then  k[1]*[1,0]/bb else k[1]*[1,0];
            D := [0];
            dZero := dZero and true;
          end if;
       else
              // State space system in form (1) with 2 first order denominator polynomials
          bb := den[1, 1] + den[2, 1] + den[1, 1]*den[2, 1] + 1;
          A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
          B[2, 1] := if abs(bb)>eps then bb else 1;
          C := if abs(bb)>eps then  k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/bb else k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])];
          D := [k[1]];
          dZero := false;
        end if;

        dss.A[1:2, 1:2] := A;
        dss.B[1:2, 1] := vector(B);
        dss.C[1, 1:2] := vector(C);
        dss.D := D;

      else
        aa := num[1,1]+1;
        bb := den[1,1]+1;
        a := -den[1, 1];
        if 1 <= n_num1 then
          // State space system in form (4)
          b := if abs(bb)>eps then bb else num[1,1]-den[1,1];
          c := if abs(bb)>eps then k[1]*(num[1, 1] - den[1, 1])/bb else k[1];
          d := k[1];
          dZero := false;
        else
         // State space system in form (5)
          b := if abs(bb)>eps then bb else if n_num1>0 then num[1,1]-den[1,1] else 1;
          c := if abs(bb)>eps then k[1]/bb else k[1];
          d := 0;
          dZero := dZero and true;
        end if;
        dss.A[1, 1] := a;
        dss.B[1, 1] := b;
        dss.C[1, 1] := c;
        dss.D[1, 1] := d;

      end if;
    /// for i=2 to degree(system)
      A[1, :] := {0,1};
      B[1, 1] := 0;
      for i in 2:max(n_den2, n_num2) loop
           // Construct state space systems of order 2
        if i <= n_den2 then
          aa := num[i, 2] + num[i, 1] +1;
          bb := den[i, 2] + den[i, 1] +1;
          A[2, :] := {-den[i, 2],-den[i, 1]};
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;

          if i <= n_num2 then
                 // State space system in form (1)

            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/bb else k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif 2*(i - n_num2) <= n_num1 then

            // State space system in form (1) with 2 first order numerator polynomials
            aa := num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] + num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] + 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]]/bb else k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif 2*(i-n_num2) -1== n_num1 then
          // State space system in form (2) with 1 first order numerator polynomial
            aa := num[2*i-n_num2-1, 1]+1;
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[2*i-n_num2-1, 1],1]/bb else k[i]*[num[2*i-n_num2-1, 1],1];
            D := [0];
            dZero := dZero and true;
          else
            // State space system in form (3)
            B[2, 1] := if abs(bb)>eps then bb else 1;
            C := if abs(bb)>eps then k[i]*[1,0]/bb else k[i]*[1,0];
            D := [0];
            dZero := dZero and true;
          end if;

        else
              // State space system in form (1) with 2 first order denominator polynomials
          bb := den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1] + den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1] + 1;
          A[2, :] := {-(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),-(den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])};
          B[2, 1] := if abs(bb)>eps then bb else 1;
  //        C := if abs(bb)>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/den[max(2*(i-n_den2)-1,i), 1]/bb else k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])];
          C := if abs(bb)>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/bb else k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])];
          D := [k[i]];
          dZero := false;
        end if;
        dss.A[2*i, 1:2*i - 2] := B[2, 1]*dss.C[1, 1:2*i - 2];
        dss.A[2*i - 1:2*i, 2*i - 1:2*i] := A;
        dss.B[2*i, 1] := if dZero then 0 else B[2, 1]*dss.D[1, 1];
        dss.C[1, 1:2*i - 2] := if dZero then fill(0, 2*i - 2) else D[1, 1]*dss.C[
          1, 1:2*i - 2];
        dss.C[1, 2*i - 1:2*i] := vector(C);
        dss.D := D*dss.D;
      end for;
   //  for i in max(2,i_d):n_den1 loop
      for i in ili:n_den1 loop
           // Construct state space systems of order 1
        bb :=  den[n_den2 + i, 1] + 1;
        a := if abs(den[n_den2 + i, 1])>eps then -den[n_den2 + i, 1] else 0.0;

        if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
              // State space system in form (4)

          c := if abs(bb)>eps then k[i_k + i]*(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/bb else 1.0;
  //        b := if abs(bb)>eps then bb else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1];
          b := if abs(bb)>eps then bb else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1]-  den[n_den2 + i, 1];
          d := k[i_k + i];
          dZero := false;
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
        // State space system in form (4)
          c := if abs(bb)>eps then k[i_k + i]*(num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/bb else 1.0;
  //        b := if abs(bb)>eps then bb else num[max(1, n_num2 + i - i_d + 1), 1];
          b := if abs(bb)>eps then bb else num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1];
          d := k[i_k + i];
          dZero := false;
        else
         // State space system in form (5)
          c := if abs(bb)>eps then k[i_k + i]/bb else k[i_k + i];
         b := if abs(bb)>eps then bb else 1;
  //        b := if abs(bb)>eps then den[n_den2 + i, 1] else 1;
          d := 0;
          dZero := dZero and true;
        end if;

  //###############################
  //        c := if abs(bb)>eps then k[1]*(num[1, 1] - den[1, 1])/bb else k[1];
  //        b := if abs(bb)>eps then bb else num[1,1]-den[1,1];
  //        d := k[1];
  //        dZero := false;
  //      else
  //       // State space system in form (5)
  //        c := if abs(bb)>eps then k[1]/bb else k[1];
  //        b := if abs(bb)>eps then bb else if n_num1>0 then num[1,1]-den[1,1] else 1;
  //        d := 0;
  //        dZero := dZero and true;
  //###############################

        dss.A[2*n_den2 + i, 1:2*n_den2 + i - 1] := b*dss.C[1, 1:2*n_den2 + i - 1];
        dss.A[2*n_den2 + i, 2*n_den2 + i] := a;
        dss.B[2*n_den2 + i, 1] := if dZero then 0 else b*dss.D[1, 1];
        dss.C[1, 1:2*n_den2 + i - 1] := if dZero then fill(0, 2*n_den2 + i - 1) else
                d*dss.C[1, 1:2*n_den2 + i - 1];
        dss.C[1, 2*n_den2 + i] := c;
        dss.D := if dZero then [0] else d*dss.D;

        end for;

      dss.C := k_total*dss.C;
      dss.D := k_total*dss.D;
    else
      dss := DiscreteStateSpace(dzp.k);
    end if;

    dss.Ts := dzp.Ts;
    dss.method := dzp.method;

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = DiscreteZerosAndPoles.Conversion<b>toDiscreteStateSpace</b>(dzp)
</pre></blockquote>

<h4>Description</h4>
<p>
This function transforms a discrete zeros-poles-gain system representation into the corresponding discrete state space representation.
To achieve well numerical condition the DiscreteZerosAndPoles transfer function is transformed into discrete state space
form by creating first and second order blocks that are connected
together in series. Every block is represented in controller
canonical form and scaled such that the gain from the input
of this block to its output is one (i.e. y(p=0) = u(p=0)),
if this is possible.
</p>
<p>
The construction of the state space blocks is the same as of the corresponding continuous transformation ( <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace\">ZerosAndPoles.Conversion.toStateSpace()</a>),
except for that the first order and second order systems are scaled by
</p>
<blockquote><pre>
       1 + n1
k_2 = --------
       1 + d1
</pre></blockquote>
<p>
and
</p>
<blockquote><pre>
       1 + n2 + n2
k_2 = -------------
       1 + d2 + d2
</pre></blockquote>
<p>
respectiively.
</p>
<p>
See <a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace\">ZerosAndPoles.Conversion.toStateSpace()</a> for more information.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q+1)/(q^2 + q + 1);

<b>algorithm</b>
  dss := Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteStateSpace(dzp);
// dss.A = [0, 1; -1, -1],
// dss.B = [0; 3],
// dss.C = [0.333, 0.333],
// dss.D = [0],
// dss.B2 = [0; 0],
</pre></blockquote>
</html>"));
  end toDiscreteStateSpace;

  end Conversion;

  encapsulated package Import
    "Package of functions to generate a DiscreteZerosAndPoles data record from imported data"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function fromFile
      "Generate a DiscreteZerosAndPoles data record by reading the polynomial coefficients or zeros and poles from a file"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DataDir;

      input String fileName=DataDir + "dzp.mat"
        "Name of the discrete zeros and poles data file"        annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
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
      output DiscreteZerosAndPoles dzp(n1=fill(0, n1), n2=fill(0, n2, 2), d1=fill(0, d1), d2=fill(0, d2, 2));
    algorithm
    //Whenever this function becomes operational the code must be rewritten if fromFile_pc2 and fromFile_zp2 are in the 'constructor'

      dzp := if ZerosAndPoles.Internal.checkRepresentation(fileName) then DiscreteZerosAndPoles.Internal.fromFile_zp( fileName) else DiscreteZerosAndPoles.Internal.fromFile_pc(
        fileName);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteZerosAndPoles.Import.fromFile(fileName)
</pre></blockquote>

<h4>Description</h4>
<p>
Reads and loads a discrete zeros-and-poles transfer function from a mat-file <code>fileName</code>.
The file must contain the sample time Ts and either the set of variables n1, n2, d1, d2, and k with
the associated first and second order polynomials or the variables p, z, and k with the poles and
zeros, written in two column arrays with real and imaginary in the first and second column respectively.
The variable k is the real gain in both cases.
</p>

<h4>Example</h4>
<blockquote><pre>
<b>algorithm</b>
  dzp := DiscreteZerosAndPoles.Import.fromFile(DataDir + &quot;/dzp.mat&quot;);
//  zp = (q^2 + 2*q + 3)/(q + 2)/(q^2 + 2*q + 2)
</pre></blockquote>
</html>"));
    end fromFile;

  function fromModel
    "Generate a DiscreteZerosAndPoles data record from a state space representation resulted from linearization of a model"

    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input String modelName "Name of the Modelica model"  annotation(Dialog(__Dymola_translatedModel(translate=true)));
    input Real T_linearize = 0
      "Point in time of simulation to linearize the model";
    input String fileName = "dslin" "Name of the result file";
    input Modelica.SIunits.Time Ts = 1 "Sample time";
    input Modelica_LinearSystems2.Utilities.Types.Method method = Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

    protected
    String fileName2 = fileName + ".mat";
    Boolean OK1 = simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
    Boolean OK2 = importInitial("dsfinal.txt");
    Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize + 1);
    Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension(fileName2, "ABCD");
    Integer nx = xuy[1];
    Integer nu = xuy[2];
    Integer ny = xuy[3];
    Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName2, "ABCD", nx + ny, nx + nu);
    String xuyName[nx + nu + ny] = readStringMatrix(fileName2, "xuyName", nx + nu + ny);

    StateSpace ss(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "Model linearized at initial point";
    DiscreteStateSpace dss(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu],
      redeclare Real B2[nx,nu]) "Model linearized at initial point";
    DiscreteStateSpace dss_siso(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,1],
      redeclare Real C[1,nx],
      redeclare Real D[1,1],
      redeclare Real B2[nx,1]) "Model linearized at initial point";

    DiscreteZerosAndPoles dummy;

    public
    output DiscreteZerosAndPoles dzp[:,:];//=fill(dummy,ny,nu);
  algorithm
    ss.A := ABCD[1:nx, 1:nx];
    ss.B := ABCD[1:nx, nx + 1:nx + nu];
    ss.C := ABCD[nx + 1:nx + ny, 1:nx];
    ss.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    ss.uNames := xuyName[nx + 1:nx + nu];
    ss.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    ss.xNames := xuyName[1:nx];

    dss := DiscreteStateSpace(ss, Ts=Ts, method=method);

    dzp := DiscreteStateSpace.Conversion.toDiscreteZerosAndPolesMIMO(dss);

  //   for ic in 1:ny loop
  //     for ib in 1:nu loop
  //       dss_siso := DiscreteStateSpace(
  //         A=dss.A,
  //         B=matrix(dss.B[:, ib]),
  //         C=transpose(matrix(dss.C[ic, :])),
  //         D=matrix(dss.D[ic, ib]),
  //         B2=matrix(dss.B2[:, ib]),
  //         Ts=dss.Ts,
  //         method=dss.method);
  //       dzp[ic, ib] := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(
  //         dss_siso);
  //     end for;
  //   end for;

  //
  // //    zp := StateSpace.Conversion.toZerosAndPolesMIMO(result);
  //     for i in 1:ny loop
  //       for j in 1:nu loop
  //         dzp[i,j] := DiscreteZerosAndPoles(zp=zp[i,j], Ts=Ts, method=method);
  //       end for;
  //     end for;

   annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteZerosAndPoles.Import.<b>fromModel</b>(
  modelName, T_linearize, fileName, Ts, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate a matrix of DiscreteZerosAndPoles data records by
linearization of a model defined by modelName. The linearization
is performed at time T_linearize of the simulation. The system is
genrated by using <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Import.fromFile\">StateSpace.Import.fromModel</a>
followed by a conversion from state space to discrete zeros-and-poles
transfer function representation.
</p>

<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
  Real T_linearize = 5;
  Modelica.SIunits.Time Ts=0.01;

<b>algorithm</b>
  dzp = Modelica_LinearSystems2.DiscreteZerosAndPoles.Import.fromModel(modelName=modelName, T_linearize=T_linearize, Ts=Ts);

//  dzp [4.48586e-006*(q^2 - 2.03695*q + 1.0381)*(q^2 - 1.9622*q + 0.963301)*(q^2 + 2*q + 1)/( (q - 1)^2*(q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal;
 //
 //    0.000897172*(q + 1)*(q^2 - 2.03695*q + 1.0381)*(q^2 - 1.9622*q + 0.963301)/( (q - 1)*(q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal;
 //
 //    -5.60444e-006*(q^2 - 1.99889*q + 1)*(q^2 + 2*q + 1)/( (q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal;
 //
 //    -0.00112089*(q - 1)*(q + 1)*(q^2 - 1.99889*q + 1)/( (q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal;
 //
 //    4.16386e-006*(q^2 - 1.99989*q + 1)*(q^2 + 2*q + 1)/( (q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal;
 //
 //    0.000832773*(q - 1)*(q + 1)*(q^2 - 1.99989*q + 1)/( (q^2 - 2.03237*q + 1.03391)*(q^2 - 1.96572*q + 0.967206) )
 //    Ts = 0.01
 //    method = Trapezoidal]
</pre></blockquote>
</html>"));
  end fromModel;

  end Import;

  encapsulated package Internal
    "Package of internal material of record DiscreteZerosAndPoles (for advanced users only)"
    extends Modelica.Icons.InternalPackage;

    import Modelica;
    import Modelica_LinearSystems2;

    function numberOfRealZeros2 "Calculate number of real zeros"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "DiscreteTransferFunction";
    output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(dtf.n)));
    algorithm
    end numberOfRealZeros2;

    function numberOfRealPoles "Calculate number of real poles"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "TransferFunction";
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
    output Real k
        "= (d+1)/(n+1), if d+1, n+1 are not zero, otherwise special cases";
  algorithm
      k := if abs(d+1) > small  and abs(n+1) > small then abs(d+1)/abs(n+1) else 1;
  end scaleFactor1;

    function scaleFactor2 "Return scale factor for second order block"
      import Modelica;
    input Real n1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real n2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real small=100*Modelica.Constants.eps;
    output Real k
        "= (d2+d1+1)/(n2+n1+1), if numerator and denominator are not zero, otherwise special cases";
    algorithm
      k := if abs(d2+d1+1) > small and abs(n2+n1+1) > small then (d2+d1+1)/(n2+n1+1) else 1;
    end scaleFactor2;

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

  encapsulated function fromFile_pc
    "Generate a DiscreteZerosAndPoles data record by reading the polynomial coefficients from a file (default file name is pc.mat)"
    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Complex;

    input String fileName="pc.mat" "Name of the zeros and poles data file"
      annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
        caption="state space system data file")));

    protected
    Integer n1n2d1d2[4]=ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(fileName);
    Integer n1=n1n2d1d2[1];
    Integer n2=n1n2d1d2[2];
    Integer d1=n1n2d1d2[3];
    Integer d2=n1n2d1d2[4];
    Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
    Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
    output DiscreteZerosAndPoles dzp(
      n1=fill(0, n1),
      n2=fill(0, n2, 2),
      d1=fill(0, d1),
      d2=fill(0, d2, 2));

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
    Real Ts[1,1] = Streams.readRealMatrix(fileName, "Ts", 1, 1);

  algorithm
    dzp.k := k;
    dzp.n1 := if n1 > 0 then n1Vector else fill(0, 0);
    dzp.n2 := if n2 > 0 then n2Matrix else fill(0, 0, 2);
    dzp.d1 := if d1 > 0 then d1Vector else fill(0, 0);
    dzp.d2 := if d2 > 0 then d2Matrix else fill(0, 0, 2);
    dzp.Ts := scalar(Ts);

  end fromFile_pc;

  encapsulated function fromFile_zp
    "Generate a DiscreteZerosAndPoles data record by reading poles and zeros from a file (default file name is zp.mat)"

    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Complex;

    input String fileName = "dzp.mat" "Name of the zeros and poles data file"
      annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
        caption="state space system data file")));
    protected
    Integer n1n2d1d2[4]=
      ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(fileName);
    Integer n1 = n1n2d1d2[1];
    Integer n2 = n1n2d1d2[2];
    Integer d1 = n1n2d1d2[3];
    Integer d2 = n1n2d1d2[4];
    Integer zSize = n1n2d1d2[1] + 2*n1n2d1d2[2];
    Integer pSize = n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
    output DiscreteZerosAndPoles dzp(
      n1=fill(0, n1),
      n2=fill(0, n2,2),
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
    Real Ts[1,1] = Streams.readRealMatrix(fileName, "Ts", 1, 1);

    Complex zeros[:] = if zSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        zerosMatrix[:, 1], zerosMatrix[:, z_2]) else fill(Complex(0), 0);
    Complex poles[:] = if pSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        polesMatrix[:, 1], polesMatrix[:, p_2]) else fill(Complex(0), 0);

  algorithm
    dzp := DiscreteZerosAndPoles(
        k=k,
        z=zeros,
        p=poles,
        Ts=scalar(Ts));
  end fromFile_zp;

    function isReal "Check of dzp whether it is a real number"
      extends Modelica.Icons.Function;

      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      input DiscreteZerosAndPoles dzp;
      output Boolean isReal;

    algorithm
      isReal := size(dzp.n1,1)+size(dzp.n2,1)+size(dzp.d1,1)+size(dzp.d2,1)==0;

    end isReal;
  end Internal;

  annotation (
    defaultComponentName="filter",
    Documentation(info="<html>
<p>
This record defines a discrete transfer function by its zeros, poles and a gain:
</p>
<pre>         product(q - z[i])
  y = k*------------------- * u
         product(q - n[i])
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

</p>
<p>
A DiscreteZeroAndPole transfer function is internally stored by the coefficients
of first and second order polynomials, and by an additional
multiplicative factor k:
</p>
<pre>         product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
  y = k*---------------------------------------------------------
         product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
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
<pre>                          (q+0.5)
  dzp = 4* -------------------------------------
           (p - 0.5)*(p - (0.6+j*0.3))*(p - (0.6-j*0.3))
</pre>
<p>
with j = Complex.j(); is defined as
</p>
<pre>
   <b>import</b> Modelica_LinearSystems2.Math.Complex;
   <b>import</b> Modelica_LinearSystems2.DiscreteZerosAndPoles;
   j = Complex.j();

   dzp = ZerosAndPoles(z = {Complex(-0.5)},
                       p = {Complex(0.5),
                            0.6+j*0.3,
                            0.6-j*0.3},
                            k=4);
</pre>
</html>"),
    Icon(
      graphics={
        Rectangle(
          lineColor={160,160,164},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid,
          extent={{-100,-100},{100,100}},
          radius=25.0),
        Text(
          lineColor={255,255,170},
          extent={{-90,-50},{90,50}},
          textString="zp")}));
end DiscreteZerosAndPoles;
