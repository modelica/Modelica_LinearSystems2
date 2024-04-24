within Modelica_LinearSystems2.Math;
operator record Polynomial "Record defining the data for a polynomial"

  Real c[:] "Polynomial coefficients (c[1]*x^n + ... c[n]*x + c[n+1])"
    annotation (Dialog);

  encapsulated package Examples
    "Package of examples to demonstrate the usage of polynomials"
    extends Modelica.Icons.ExamplesPackage;
    import Modelica;

    function polynomialOperations1
      "Demonstrate basic algebraic operations on polynomials"
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2.Math.Polynomial;

      output Boolean ok;

    protected
      Polynomial t=Polynomial.x();
      Polynomial p;
    algorithm
      p := t^2 + 8*t + 16;

      print("y1 = " + String(p, name="t"));
      print("y2 = " + String(p*p, name="t"));
      ok := true;
    end polynomialOperations1;

    function polynomialOperations2
      "Demonstrate differentiation and integration of polynomials"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica.Utilities.Streams.print;

      output Boolean ok;
    protected
      Polynomial p=Polynomial({6,4,-3}) "6*x^2 + 4*x - 3";
      Polynomial int_p;
      Polynomial der_p;
      TransferFunction tf;
      Real x;
      Real int_val1;
      Real int_val2;
      Real der_val1;
      Real der_val2;
      Real der_val3;
      Real der_val4;
    algorithm
      print("Show integration and differentation of y = " + String(p));
      int_p := Polynomial.integral(p);
      der_p := Polynomial.derivative(p);

      print("  Integral of polynomial  : " + String(int_p));
      print("  Derivative of polynomial: " + String(der_p));

      x := 2;
      print(
        "Compute derivatives and integral directly and via polynomials" +
        " above (for p(x=" + String(x) + ")):");
      int_val1 := Polynomial.integralValue(p, x);
      int_val2 := Polynomial.evaluate(int_p, x);
      print("      integralValue(p,x) = " + String(int_val1) +
        ", evaluate(integral(p),x) = " + String(int_val2));

      der_val1 := Polynomial.derivativeValue(p, x);
      der_val2 := Polynomial.evaluate(der_p, x);
      print("    derivativeValue(p,x) = " + String(der_val1) +
        ", evaluate(derivative(p),x) = " + String(der_val2));
      der_val3 := Polynomial.derivativeValue(p, x, 2);
      der_val4 := Polynomial.evaluate(Polynomial.derivative(der_p), x);
      print("  derivativeValue(p,x,2) = " + String(der_val3) +
        ", evaluate(derivative(derivative(p)),x) = " + String(der_val4));
      tf := der_p/int_p;
      print("  derivative/integral of polynomial = " + String(tf, name="x"));
      ok := true;
    end polynomialOperations2;

    function PascalTriangle "Generate and print Pascals triangle"
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica.Utilities.Streams.print;

      output Boolean ok;
    protected
      Polynomial p=Polynomial.x();
      Polynomial p1;
    algorithm
      p1 := p + 1;
      print("Pascals triangle up to order 20:");
      for i in 0:20 loop
        print("  " + String(p1^i));
      end for;
      ok := true;
    end PascalTriangle;

    function plotPolynomial
      "Demonstrate plotting of polynomial with automatic range and legend determination"

      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Complex;

      output Boolean ok;

    protected
      Complex c[:]={-2 + 0*j,Complex(0),2 + 0*j,7 + j,7 - j};
      Polynomial p=Polynomial(c);
    algorithm
      Polynomial.plot(p);
      ok := true;
      annotation (__Dymola_interactive=true);
    end plotPolynomial;

    model MoveMassByPolynomial "Demonstrate meaning of polynom's derivatives"
      extends Modelica.Icons.Example;

      Utilities.MoveByPolynomial positionByPolynomial(useSupport=false, c3=30) annotation (Placement(transformation(extent={{-70,10},{-50,30}})));
      Modelica.Mechanics.Translational.Components.Mass mass1(m=1) annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
      Modelica.Mechanics.Translational.Components.Mass mass2(
        m=1,
        stateSelect=StateSelect.always,
        s(fixed=true, start=0.1),
        v(fixed=true)) annotation (Placement(transformation(extent={{20,10},{40,30}})));
      Modelica.Mechanics.Translational.Components.SpringDamper springDamper(c=5, d=0.1) annotation (Placement(transformation(extent={{-10,10},{10,30}})));
      Modelica.Mechanics.Translational.Sensors.PositionSensor positionSensor2 annotation (Placement(transformation(extent={{50,10},{70,30}})));
      Utilities.EvalPolynomial3 evalPolynomial3(coefficients={positionByPolynomial.c3,0,0,0}) annotation (Placement(transformation(extent={{-60,-90},{-40,-70}})));
      Modelica.Mechanics.Translational.Sensors.PositionSensor positionSensor1 annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
      Modelica.Mechanics.Translational.Sensors.AccSensor accSensor annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
      Modelica.Blocks.Math.Add add[3](each k2=-1) annotation (Placement(transformation(extent={{40,-80},{60,-60}})));
    equation
      connect(mass1.flange_a, positionByPolynomial.flange) annotation (Line(points={{-40,20},{-50,20}},
                                                                                                      color={0,127,0}));
      connect(mass1.flange_b, springDamper.flange_a) annotation (Line(points={{-20,20},{-10,20}},
                                                                                                color={0,127,0}));
      connect(springDamper.flange_b,mass2. flange_a) annotation (Line(points={{10,20},{20,20}},
                                                                                              color={0,127,0}));
      connect(mass2.flange_b, positionSensor2.flange) annotation (Line(points={{40,20},{50,20}}, color={0,127,0}));
      connect(positionSensor2.s, positionByPolynomial.u) annotation (Line(points={{71,20},{80,20},{80,-10},{-80,-10},{-80,20},{-72,20}}, color={0,0,127}));
      connect(positionSensor2.s, evalPolynomial3.u) annotation (Line(points={{71,20},{80,20},{80,-10},{-80,-10},{-80,-80},{-62,-80}}, color={0,0,127}));
      connect(mass1.flange_b, positionSensor1.flange) annotation (Line(points={{-20,20},{-16,20},{-16,-20},{-10,-20}}, color={0,127,0}));
      connect(mass1.flange_b, speedSensor.flange) annotation (Line(points={{-20,20},{-16,20},{-16,-40},{-10,-40}}, color={0,127,0}));
      connect(mass1.flange_b, accSensor.flange) annotation (Line(points={{-20,20},{-16,20},{-16,-60},{-10,-60}}, color={0,127,0}));
      connect(evalPolynomial3.y, add.u2) annotation (Line(points={{-39,-80},{20,-80},{20,-76},{38,-76}}, color={0,0,127}));
      connect(positionSensor1.s, add[1].u1) annotation (Line(points={{11,-20},{22,-20},{22,-64},{38,-64}}, color={0,0,127}));
      connect(speedSensor.v, add[2].u1) annotation (Line(points={{11,-40},{20,-40},{20,-64},{38,-64}}, color={0,0,127}));
      connect(accSensor.a, add[3].u1) annotation (Line(points={{11,-60},{18,-60},{18,-64},{38,-64}}, color={0,0,127}));
      annotation (
        experiment(StopTime=10),
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false), graphics={Text(
              extent={{-80,64},{80,44}},
              textColor={28,108,200},
              textString="Note: positionByPolynomial.c3 must be less then 1/mass2.s.start^2"),
            Line(
              points={{-40,50},{-54,32}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Open}),
            Line(
              points={{58,50},{34,30}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Open})}),
        Documentation(info="<html>
<p>
The motion of the <code>mass1</code> is forced depending on the motion
of <code>mass2</code>. The dependency is given by a&nbsp;cubic polynomial.
</p>
<p>
This example demonstrates the meaning of the function evaluating
the 2nd derivative of a&nbsp;polynomial with respect to time.
Without this function, the DAE index reduction of the mechanical
system would fail since <code>positionByPolynomial.s</code>, which prescribes
the motion of <code>mass1</code>, needs to be differentiated twice.
</p>
<p>
In order to check proper differentiation, the block <code>evalPolynomial3</code>
is provided which only evaluates a&nbsp;cubic polynomial of the same
coefficients like <code>positionByPolynomial</code>.
The signal differences given by <code>add</code> block shall, thus, be zero
(see simulation results of <code>add[1].y</code>, <code>add[2].y</code>
and <code>add[3].y</code>).
</p>
</html>"));
    end MoveMassByPolynomial;

    package Utilities "Utility classes for polynomial examples"
      extends Modelica.Icons.UtilitiesPackage;

      model MoveByPolynomial "Forced movement of a flange according to a given cubic polynomial"
        //extends Modelica.Mechanics.Translational.Interfaces.PartialCompliantWithRelativeStates;
        extends Modelica.Mechanics.Translational.Interfaces.PartialElementaryOneFlangeAndSupport2;
        import Modelica_LinearSystems2.Math.Polynomial;
        import Modelica.Units.SI;

        parameter Real c3(start=1) "Coefficient c3 of cubic polynomial";

        Modelica.Blocks.Interfaces.RealInput u "Polynomial's variable" annotation (Placement(transformation(
              extent={{-20,-20},{20,20}},
              origin={-120,0})));

      protected
        Polynomial p1 = Polynomial({c3,0,0,0}) "Cubic polynomial";

      equation
        s = Polynomial.evaluate(p1, u);

        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false),
            graphics={
              Rectangle(
                extent={{0,20},{100,-20}},
                lineColor={0,127,0},
                fillColor={160,215,160},
                fillPattern=FillPattern.Solid),
              Line(points={{11,32},{70,32}}, color={0,127,0}),
              Line(points={{40,52},{40,32}}, color={0,127,0}),
              Line(points={{10,-32},{70,-32}}, color={0,127,0}),
              Line(points={{40,-72},{40,-80},{0,-80},{0,-100}}, color={0,127,0}),
              Line(points={{-50,60},{-50,-60}}, color={192,192,192}),
              Line(points={{-90,0},{-10,0}}, color={192,192,192}),
              Line(
                points={{-88,-52},{-84,-30},{-72,-6},{-50,0},{-28,6},{-16,30},{-12,52}},
                color={0,0,127},
                smooth=Smooth.Bezier),
              Text(
                extent={{-150,100},{150,60}},
                textColor={0,0,255},
                textString="%name"),
              Text(
                extent={{-150,-41},{150,-71}},
                textColor={0,0,0},
                textString="c3=%c3"),
              Line(points={{40,-32},{40,-38}},
                color={0,127,0})}),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<p>
The <code>flange</code> is <strong>forced</strong> to move
relative to the <code>support</code> connector with a&nbsp;predefined motion
according to the input signal&nbsp;<code>u</code> as follows.
</p>
<blockquote><pre>
position of flange:     s = c3 * u^3
velocity of flange:     v = c3 * 3 * u^2
acceleration of flange: a = c3 * 6 * u
</pre></blockquote>
<p>
whereby the evaluation is done by the record
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial\">Polynomial</a>.
</p>
</html>"));
      end MoveByPolynomial;

      model EvalPolynomial3 "Evaluate cubic polynomial in time domain"
        extends Modelica.Blocks.Interfaces.SIMO(final nout=3);
        import Modelica_LinearSystems2.Math.Polynomial;

        parameter Real coefficients[4] = {1,0,0,0} "Coefficients of cubic polynomial";
      protected
        Polynomial p1 = Polynomial(coefficients) "Cubic polynomial";
        Polynomial der_p1 = Polynomial.derivative(p1);
        Real du "Derivative 1st of u";

      equation
        du = der(u);
        y[1] = Polynomial.evaluate(p1, u);
        y[2] = Polynomial.evaluate(der_p1, u) * du;
        y[3] = Polynomial.evaluate(Polynomial.derivative(der_p1), u) * du * du
             + Polynomial.evaluate(der_p1, u) * der(du);
        annotation (Icon(graphics={
              Line(points={{0,80},{0,-80}},     color={192,192,192}),
              Line(points={{-80,0},{80,0}},  color={192,192,192}),
              Line(
                points={{-60,-80},{-50,-20},{-30,4},{0,0},{30,-4},{50,20},{60,80}},
                color={0,0,127},
                smooth=Smooth.Bezier)}), Documentation(info="<html>
<p>
Evaluate a&nbsp;cubic polynomial <var>p</var>(&nbsp;<var>u</var>(<var>t</var>)&nbsp;)
and its derivatives in time domain. The output
vector&nbsp;<code>y</code> is composed as follows.
</p>
<blockquote><pre>
y[1] = c[1] * u(t)^3 + c[2] * u(t)^2 + c[3] * u(t) + c[4],
y[2] = der(y[1]),
y[3] = der(y[2]),
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
c[i] := coefficients[i].
</pre></blockquote>
</html>"));
      end EvalPolynomial3;
      annotation (Documentation(info="<html>
<p>
This library contains different utility components
for examples involving Polynomial record. Usually, there is no need to
use models in this library directly.
</p>
</html>"));
    end Utilities;
  end Examples;

  encapsulated operator 'constructor'
    "Collection of operators to construct a Polynomial data record"
    import Modelica;

    function fromVector
      "Generate a Polynomial data record from a vector of coefficients"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Real c[:] "Coefficients of polynomial in descending order";
      output Polynomial p(redeclare Real c[size(c, 1)]) "Polynomial";
    algorithm
      p.c := c;
      annotation (Documentation(info="<html>
<p>
This function constructs a&nbsp;polynomial from given coefficients.
</p>

<h4>
Example:
</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

Polynomial( {1/4, 8, 0, 0, 2} );
// 0.25*x^4 + 8*x^3 + 2
</pre></blockquote>
</html>"));
    end fromVector;

    function fromReal "Generate a Polynomial data record from a real value"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Real r "Value of Real variable";
      output Polynomial p(redeclare Real c[1]);
    algorithm
      p.c := {r};
      annotation (Documentation(info="<html>
<p>
This function constructs a&nbsp;polynomial from a&nbsp;real value.
</p>

<h4>
Example:
</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

p1 = Polynomial(4);
p1;
// 4
</pre></blockquote>
</html>"));
    end fromReal;

    function fromZeros "Generate a Polynomial data record from given zeros"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica.Utilities.Streams;
      import Complex;

      input Complex roots[:]
        "Zeros of polynomial (must be real or conjugate complex pairs)";
      output Polynomial p(redeclare Real c[size(roots, 1) + 1])
        "Polynomial that corresponds to the zeros";

    protected
      Integer nr=size(roots, 1) "Number of roots";
      Integer nc=nr + 1 "Length of coefficient vector";
      Integer nn;
      Integer i;
    algorithm
      p.c[nc] := 1;
      i := 1;
      nn := 1;
      while i <= nr loop
        // Determine whether zero is real or conjugate complex
        if roots[i].im == 0.0 then
          // real root
          p.c := Polynomial.Internal.mult(
                p.c,
                nn,
                {1,-roots[i].re},
                nc);

          i := i + 1;
          nn := nn + 1;
        else
          // complex root; check that it is a conjugate complex pair
          assert(i < nr, "Roots do not define a real valued polynomial\n" +
            "(roots[" + String(nr) + "] is complex without complex conjugate)");
          assert(roots[i].re == roots[i + 1].re,
            "No conjugate complex zero pair\n" + "  roots[" + String(i) +
            "].re = " + String(roots[i].re) + "\n" + "  roots[" + String(i + 1)
             + "].re = " + String(roots[i].re) + "\n" +
            "and the two values should be identical since conjugate complex pair required.");
          assert(roots[i].im == -roots[i + 1].im,
            "No conjugate complex zero pair\n" + "  roots[" + String(i) +
            "].im = " + String(roots[i].im) + "\n" + "  roots[" + String(i + 1)
             + "].im = " + String(roots[i + 1].im) + "\n" +
            "and the two values should be identical with opposite sign, since conjugate complex pair required.");
          p.c := Polynomial.Internal.mult(
                p.c,
                nn,
                {1,-2*roots[i].re,roots[i].re^2 + roots[i].im^2},
                nc);
          i := i + 2;
          nn := nn + 2;
        end if;
      end while;
      annotation (Documentation(info="<html>
<p>
This function constructs a&nbsp;polynomial from given zeros
(also called roots). The zeros are defined as a&nbsp;vector
of Complex numbers. Since only polynomials with real coefficients are supported,
complex zeros must be defined as conjugate complex pairs.
It is required that complex conjugate pairs must directly
follow each other. An error occurs if this is not the case.
</p>

<h4>
Example:
</h4>
<p>
The polynomial
</p>
<blockquote>
<var>y</var> = (<var>s</var> - 1) * ( <var>s</var> - (2+3j) ) * (<var>s</var> - (2-3j))
</blockquote>
<p>
with j=sqrt(-1), is defined as
</p>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial
Polynomial( {Complex(1),
             Complex(2,  3),
             Complex(2, -3)} );
// x^3 - 5*x^2 + 17*x - 13
</pre></blockquote>
</html>"));
    end fromZeros;

    annotation (Documentation(info="<html>
<p>
This package contains the default constructors for a&nbsp;data
record of polynomial. Usually, the functions shall not be accessed
directly.
</p>
</html>"));
  end 'constructor';

  encapsulated operator '-'
    "Collection of operators for subtraction of polynomials"
    import Modelica;

    encapsulated function negate "Unary minus (multiply polynomial by -1)"

      import Modelica_LinearSystems2.Math.Polynomial;
      input Polynomial p;
      output Polynomial result(redeclare Real c[size(p.c, 1)]) "= -p";
    algorithm
      result.c := -p.c;
    end negate;

    encapsulated function subtract "Subtract two polynomials (p1 - p2)"
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial p1;
      input Polynomial p2;
      output Polynomial result(redeclare Real c[max(size(p1.c, 1), size(p2.c, 1))])
        "= p1 - p2";
    algorithm
      // Auxiliary variables not used, to enforce function inlining
      result.c := cat(
            1,
            zeros(max(size(p1.c, 1), size(p2.c, 1)) - size(p1.c, 1)),
            p1.c) - cat(
            1,
            zeros(max(size(p1.c, 1), size(p2.c, 1)) - size(p2.c, 1)),
            p2.c);
    end subtract;
    annotation (Documentation(info="<html>
<p>
This package contains operators for subtraction of Polynomial data records.
</p>
</html>"));
  end '-';

  encapsulated operator function '+' "Add two polynomials (p1 + p2)"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p1;
    input Polynomial p2;
    output Polynomial result(redeclare Real c[max(size(p1.c, 1), size(p2.c, 1))])
      "= p1 + p2";
  algorithm
    // Auxiliary variables not used, to enforce function inlining
    result.c := cat(
        1,
        zeros(max(size(p1.c, 1), size(p2.c, 1)) - size(p1.c, 1)),
        p1.c) + cat(
        1,
        zeros(max(size(p1.c, 1), size(p2.c, 1)) - size(p2.c, 1)),
        p2.c);
  end '+';

  encapsulated operator function '*' "Multiply two polynomials (p1 * p2)"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p1;
    input Polynomial p2;
    output Polynomial result(redeclare Real c[size(p1.c, 1) + size(p2.c, 1) - 1])
      "= p1 * p2";
  protected
    Integer n1=size(p1.c, 1);
    Integer n2=size(p2.c, 1);
    Integer n3=n1 + n2 - 1;
    Real ck;
  algorithm
    for k in 1:n3 loop
      ck := 0.0;
      for j in max(1, k + 1 - n2):min(k, n1) loop
        ck := ck + p1.c[j]*p2.c[k + 1 - j];
      end for;
      result.c[k] := ck;
    end for;
  end '*';

  encapsulated operator function '/' "Divide two polynomials (p1 / p2)"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.TransferFunction;

    input Polynomial p1;
    input Polynomial p2;

    output TransferFunction tf(n=p1.c, d=p2.c);
    //only for tfpoly
  algorithm

    assert(size(p2.c, 1) > 0,
      "Denominator polynomial p2 must have at least one element, however\n" +
      "denominator is an empty polynomial. This is not allowed for p1/p2.");

  end '/';

  encapsulated operator function '^' "Integer power of polynomial (p^n)"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica.Utilities.Streams.print;

    input Polynomial p;
    input Integer n(min=0) = 1 "p^n shall be computed";
    output Polynomial result(redeclare Real c[max((size(p.c, 1) - 1)*n + 1, 1)])
      "= p^n";
  protected
    Integer n_p=size(p.c, 1);
    Integer n_power_p=max((n_p - 1)*n + 1, 1);
  algorithm
    if n == 0 then
      result.c := {1};
    else
      result.c[n_power_p - n_p + 1:n_power_p] := p.c;
      for i in 2:n loop
        result.c := Polynomial.Internal.mult(
            result.c,
            (n_p - 1)*(i - 1) + 1,
            p.c,
            n_power_p);
      end for;
    end if;
  end '^';

  encapsulated operator function '=='
    "Check whether two polynomials are identical"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p1;
    input Polynomial p2;
    input Real eps(min=0) = 0
      "Two coefficients c1 and c2 of the two polynomials are identical if abs(c1-c2) <= eps";
    output Boolean same "=true, if identical";
  protected
    Integer n1=size(p1.c, 1);
    Integer n2=size(p2.c, 1);
  algorithm
    if n1 == n2 then
      same := true;
      for i in 1:n1 loop
        if abs(p1.c[i] - p2.c[i]) > eps then
          same := false;
        end if;
      end for;
    else
      same := false;
    end if;
  end '==';

  encapsulated operator function 'String'
    "Transform Polynomial into a String representation"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica;

    input Polynomial p
      "Polynomial to be transformed in a String representation";
    input Integer significantDigits=6
      "Number of significant digits that are shown";
    input String name="x" "Independent variable name used for printing";
    output String s="";
  protected
    Boolean outputCoefficient;
    Integer power;
    Real ci;
    String v;
    Integer n=size(p.c, 1);
    Boolean first=true;
  algorithm
    if n == 0 then
      s := "0";
    else
      for i in 1:n loop
        if p.c[i] <> 0 or i == n then
          power := n - i;
          ci := p.c[i];
          if first then
            first := false;
          else
            if ci > 0 then
              s := s + " + ";
            elseif ci < 0 then
              s := s + " - ";
              ci := abs(ci);
            end if;
          end if;

          outputCoefficient := power == 0 or abs(ci - 1) > Modelica.Constants.eps;
          if outputCoefficient then
            s := s + String(ci, significantDigits=significantDigits);
          end if;
          if outputCoefficient and power >= 1 then
            s := s + "*";
          end if;
          if name == "" then
            v := "?";
          else
            v := name;
          end if;
          if power >= 2 then
            s := s + v + "^" + String(power);
          elseif power == 1 then
            s := s + v;
          end if;
        end if;
      end for;
    end if;
  end 'String';

  encapsulated function x "Generate a base polynomial y=x"
    import Modelica_LinearSystems2.Math.Polynomial;

    output Polynomial p(redeclare Real c[2]) "Polynomial";
  algorithm
    p.c := {1,0};
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
p = Polynomial.<strong>x</strong>();
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a&nbsp;base polynomial.
</p>

<h4>
Example:
</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();      // Generate a base polynomial
p1 = -6*x^2 + 4*x -3;
p1;
// -6*x^2 + 4*x - 3
</pre></blockquote>
</html>"));
  end x;

  encapsulated function fitting
    "Compute a Polynomial that fits a set of data points in a least-squares sense"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;

    input Real x[:] "Abscissa data values";
    input Real y[size(x, 1)] "Ordinate data values";
    input Integer order(min=1)
      "Order of desired polynomial that fits the data points (x,y)";
    output Polynomial p(redeclare Real c[order + 1])
      "Polynomial that fits the date points in a least squares sense";

  protected
    Real V[size(x, 1), order + 1] "Vandermonde matrix";
  algorithm
    // Construct Vandermonde matrix
    V[:, order + 1] := ones(size(x, 1));
    for j in order:-1:1 loop
      V[:, j] := {x[i]*V[i, j + 1] for i in 1:size(x, 1)};
    end for;

    // Solve least squares problem
    p.c := Modelica.Math.Matrices.leastSquares(V, y);
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
p = Polynomial.<strong>fitting</strong>(x, y, order);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes a&nbsp;Polynomial <var>y</var>&nbsp;=&nbsp;<var>p</var>(<var>x</var>)
of degree &quot;order&quot; that fits the data
<var>p</var>(<var>x</var><sub><var>i</var></sub>)&nbsp;-&nbsp;<var>y</var><sub><var>i</var></sub>
for <var>i</var>&nbsp;=&nbsp;1,&nbsp;..&nbsp;<var>n</var><sub><var>x</var></sub>
in a&nbsp;least squares sense.
</p>
</html>"));
  end fitting;

  encapsulated function degree "Return degree of polynomial"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p;
    output Integer result "Degree of polynomial p";
  protected
    Integer s;
  algorithm
    s := size(p.c, 1);
    for i in 1:s loop
      // added correct code for degree calculation
      if p.c[i] <> 0 then
        result := s - i;
        break;
      end if;
    end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Polynomial.<strong>degree</strong>(p);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the order of a&nbsp;Polynomial <var>p</var>.
</p>

<h4>
Example:
</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

Polynomial p = Polynomial({0,0,2,1});
// 2*x + 1

Polynomial.degree(p)
//  = 1
</pre></blockquote>
</html>"));
  end degree;

  encapsulated function plot "Plot polynomial y=p(x)"

    import Modelica;
    import DymolaCommands.Plot;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica.Utilities.Strings;
    import Complex;

    input Polynomial p "Polynomial to be plotted";
    input Integer nPoints(min=2) = 200 "Number of points";
    input Boolean autoLabel=true "True, if automatically selected labels";
    input String xLabel="" "Abszissa description, if autoLabel = false";
    input String yLabel="" "Ordinate description, if autoLabel = false";
    input Boolean autoRange=true
      "True, if abszissa range is automatically determined";
    input Real x_min=-1.0 "Minimum abszissa value, if autoRange = false";
    input Real x_max=1.0 "Maximum abszissa value, if autoRange = false";
  protected
    Real x_min2;
    Real x_max2;
    Real dx;
    Real x[nPoints];
    Real y[nPoints];
    Boolean OK;
    String yLabel2;
    Complex p_zeros[:];
    Complex pd_zeros[:];
    Complex points[:];
    String argument="x";
  algorithm
    /* Determine suitable x_min and x_max:
        Plotted range should contain roots and extrema
     */

    if autoRange then

      x_min2 := -1;
      x_max2 := +1;
      if size(p.c, 1) > 1 then
        p_zeros := Polynomial.roots(p);
        pd_zeros := Polynomial.roots(Polynomial.derivative(p));
        points := cat(
            1,
            p_zeros[:],
            pd_zeros[:]);

        x_min2 := points[1].re;
        x_max2 := points[1].re;
        for i in 2:size(points, 1) loop
          x_min2 := min(x_min2, points[i].re);
          x_max2 := max(x_max2, points[i].re);
        end for;

        if x_max2 > x_min2 then
          x_min2 := x_min2 - 0.1*(x_max2 - x_min2);
          x_max2 := x_max2 + 0.1*(x_max2 - x_min2);
        else
          x_min2 := x_min2 - 1;
          x_max2 := x_max2 + 1;
        end if;
      end if;
    else
      x_min2 := x_min;
      x_max2 := x_max;
    end if;

    // Compute polynomial values
    dx := (x_max2 - x_min2)/(nPoints - 1);
    x := x_min2:dx:x_max2;
    y := Polynomial.evaluate(p, x);

    // Determine labels
    if autoLabel then
      yLabel2 := String(p);
      if Strings.length(yLabel2) >= 100 then
        yLabel2 := Strings.substring(yLabel2, 1, 96) + " ...";
      end if;
    else
      yLabel2 := yLabel;
    end if;

    // Plot Polynomial
    OK := Plot.plotArray(x, y, legend=yLabel2);

  equation

    annotation (
      __Dymola_interactive=true,
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Polynomial.<strong>plot</strong>(
  p, nPoints, autoLabel, xLabel, yLabel, autoRange, x_min, x_max);
</pre></blockquote>

<h4>Description</h4>
<p>
Plot the given polynomial. If default arguments are used, as in:
</p>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

p = Polynomial({1,2,3});
Polynomial.plot(p);
</pre></blockquote>
<p>
then the abszissa range is determined in such a&nbsp;way that all
roots (i.e., <var>p</var>(<var>x</var>)&nbsp;=&nbsp;0) and all
extrema (i.e, <var>dp</var>(<var>x</var>)&nbsp;= &part;<var>p</var>/&part;<var>x</var>&nbsp;=&nbsp;0)
are in the plotted range. As default legend, the String representation
of the polynomial is used as generated by operator-&apos;String&apos;
function (called by <code>String(p)</code>).
</p>
</html>"));
  end plot;

  encapsulated function derivative "Derivative of polynomial"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p;
    output Polynomial der_p(redeclare Real c[size(p.c, 1) - 1])
      "Derivative of Polynomial p";
  protected
    Integer n=size(p.c, 1);
  algorithm
    for j in 1:n - 1 loop
      der_p.c[j] := p.c[j]*(n - j);
    end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
der_p = Polynomial.<strong>derivative</strong>(p);
</pre></blockquote>

<h4>Description</h4>
<p>
Define a&nbsp;polynomial <code>der_p</code> being a&nbsp;first partial derivative of
the polynomial&nbsp;<code>p</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

der_p = Polynomial.derivative(p1);
// -12*x + 4
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.integral\">Polynomial.integral</a>,
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.derivativeValue\">Polynomial.derivativeValue</a>.
</p>
</html>"));
  end derivative;

  encapsulated function integral "Indefinite integral of polynomial"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial";
    output Polynomial integral_p(redeclare Real c[size(p.c, 1) + 1])
      "Indefinite integral of polynomial p";
  protected
    Integer n=size(p.c, 1) + 1 "Number of coefficients of integral";
  algorithm
    for j in 1:n - 1 loop
      integral_p.c[j] := p.c[j]/(n - j);
    end for;
    integral_p.c[n] := 0;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
integral_p = Polynomial.<strong>integral</strong>(p);
</pre></blockquote>

<h4>Description</h4>
<p>
Define a&nbsp;polynomial <code>integral_p</code> being an indefinite integral of
the polynomial&nbsp;<code>p</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

int_p = Polynomial.integral(p1);
// -2*x^3 + 2*x^2 - 3*x
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.derivative\">Polynomial.derivative</a>.
</p>
</html>"));
  end integral;

  encapsulated function evaluate
    "Evaluate a polynomial at a given Real abszissa value"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial to be evaluated";
    input Real x "Abszissa value";
    output Real y "Value of polynomial at x";

  protected
    Integer n=size(p.c, 1);
  algorithm
    y := p.c[1];
    for j in 2:n loop
      y := p.c[j] + x*y;
    end for;
    annotation (
      derivative(zeroDerivative=p) = Polynomial.Internal.evaluate_der,
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
y = Polynomial.<strong>evaluate</strong>(p, x);
</pre></blockquote>

<h4>Description</h4>
<p>
Evaluate a&nbsp;polynomial&nbsp;<code>p</code> at a&nbsp;given
<em>Real</em> abszissa value&nbsp;<code>x</code>.
</p>

<h4>Note</h4>
<p>
For this function, functions for the first and second derivatives
are provided, so that the function can be seamlessly used e.g. in
mechanical models.
<!-- , cf. example . -->
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.evaluate(p1, 1);
// = -5.0
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.derivativeValue\">Polynomial.derivativeValue</a>,
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.integralValue\">Polynomial.integralValue</a>,
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Polynomial.evaluateComplex</a>.
</p>
</html>"));
  end evaluate;

  encapsulated function evaluateMatrix
    "Evaluate a polynomial with a matrix argument"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial to be evaluated";
    input Real X[:, size(X, 1)] "Square matrix argument";
    output Real Y[size(X, 1), size(X, 2)] "Value of polynomial at X";

  protected
    Integer n=size(p.c, 1);
  algorithm
    Y := zeros(size(X, 1), size(X, 2));
    for j in 1:n loop
      Y := X*Y + diagonal(p.c[j]*ones(size(X, 1)));
    end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Y = Polynomial.<strong>evaluateMatrix</strong>(p, X);
</pre></blockquote>

<h4>Description</h4>
<p>
Evaluates the given polynomial <var>p</var> of order <var>n</var> with its coefficients
<var>c</var><sub><var>i</var></sub> so that
</p>
<blockquote><pre>
<strong>Y</strong> = p.c[1]*<strong>X</strong>^n + p.c[2]*<strong>X</strong>^(n-1) + ... + p.c[n]*<strong>X</strong> + p.c[n+1]*<strong>I</strong>
</pre></blockquote>

<h4>Note</h4>
<p>
The matrix <strong>X</strong> must be square.
Horner's method is used for polynomial evaluation.
</p>
</html>"));
  end evaluateMatrix;

  encapsulated function evaluateComplex
    "Evaluate a polynomial at a given Complex abszissa value"
    import Complex;
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial to be evaluated";
    input Complex x "Complex abszissa value";
    output Complex y "Complex value of polynomial at x";
  protected
    Integer n=size(p.c, 1);
  algorithm
    // Horners scheme for complex numbers
    y := Complex(p.c[1], 0);
    for j in 2:n loop
      y := Complex(p.c[j], 0) + x*y;
    end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
y = Polynomial.<strong>evaluateComplex</strong>(p, x);
</pre></blockquote>

<h4>Description</h4>
<p>
Evaluate a&nbsp;polynomial&nbsp;<code>p</code> at a&nbsp;given
<em>Complex</em> abszissa value&nbsp;<code>x</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.evaluateComplex(p1, Complex(2, 3));
// 35 - 60*j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluate\">Polynomial.evaluate</a>.
</p>
</html>"));
  end evaluateComplex;

  encapsulated function derivativeValue
    "Evaluate an i-the derivative of polynomial at a given Real abszissa value"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica.Utilities.Streams.print;

    input Polynomial p "Polynomial";
    input Real x "Abszissa value";
    input Integer i(min=1) = 1
      "i-th derivative to be evaluated; i=1 is first derivative";
    output Real der_y "Value of i-th derivative of polynomial at x";
    // annotation(derivative(zeroDerivative=p, zeroDerivative=i)=Polynomial.Internal.derivativeValue_der);
  protected
    Integer n=size(p.c, 1);
  algorithm
    if i > n - 1 then
      der_y := 0.0;
    else
      der_y := p.c[1]*product(n - k for k in 1:i);
      for j in 2:n - i loop
        der_y := p.c[j]*product(n - j + 1 - k for k in 1:i) + x*der_y;
      end for;
    end if;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
der_y = Polynomial.<strong>derivativeValue</strong>(p, x, i);
</pre></blockquote>

<h4>Description</h4>
<p>
Evaluate an <code>i</code>-th derivative of a&nbsp;polynomial&nbsp;<code>p</code>
at a&nbsp;given <em>Real</em> abszissa value&nbsp;<code>x</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.derivativeValue(p1, 1, 2);
// = -12.0
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluate\">Polynomial.evaluate</a>,
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.integralValue\">Polynomial.integralValue</a>.
</p>
</html>"));
  end derivativeValue;

  encapsulated function integralValue
    "Evaluate integral of polynomial from x_low to x_high"
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial to be integrated";
    input Real x_high "High integrand value";
    input Real x_low=0 "Low integrand value, default 0";
    output Real integral=0.0 "Integral of polynomial p from x_low to x_high";

  protected
    Integer n=size(p.c, 1) "Order of integrated polynomial";
    Real y_low=0 "value at lower integrand";
  algorithm
    for j in 1:n loop
      integral := x_high*(p.c[j]/(n - j + 1) + integral);
      y_low := x_low*(p.c[j]/(n - j + 1) + y_low);
    end for;
    integral := integral - y_low;
    annotation (
      derivative(zeroDerivative=p) = Polynomial.Internal.integralValue_der,
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
integral = Polynomial.<strong>integralValue</strong>(p, x_high, x_low);
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate value of integral of a&nbsp;polynomial&nbsp;<code>p</code>
at a&nbsp;given range between&nbsp;<code>x_low</code> and&nbsp;<code>x_high</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.integralValue(p1, 2);
// = -14.0

Polynomial.integralValue(p1, 2, -2);
// = -44.0
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluate\">Polynomial.evaluate</a>,
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.derivativeValue\">Polynomial.derivativeValue</a>.
</p>
</html>"));
  end integralValue;

  encapsulated function roots
    "Determine zeros of polynomial, i.e., points x with p(x)=0"
    import Complex;
    import Modelica_LinearSystems2.ComplexMathAdds.Vectors;
    import Modelica_LinearSystems2.Math.Polynomial;

    input Polynomial p "Polynomial";
    input Boolean printRoots=false "True, if roots shall be pretty printed";
    output Complex result[:]=fill(Complex(0, 0), Polynomial.numberOfRoots(p))
      "Roots of polynomial";

  algorithm
    result := Polynomial.rootsOfNonZeroHighestCoefficientPolynomial(p,
      Polynomial.numberOfRoots(p));
    if printRoots then
      Vectors.print("", result);
    end if;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Polynomial.<strong>roots</strong>(p, printRoots);
</pre></blockquote>

<h4>Description</h4>
<p>
The roots of the given polynomial&nbsp;<code>p</code> are determined
and are returned as a&nbsp;vector of Complex elements.
If the optional argument <code>printRoots</code>
is set to true, the roots are also printed at once.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.roots(p1);
// {0.333333 + 0.62361*j, 0.333333 - 0.62361*j}

Polynomial.roots(p1, true);
//  =
//    0.333333 + 0.62361*j
//    0.333333 - 0.62361*j
// {0.333333 + 0.62361*j, 0.333333 - 0.62361*j}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.numberOfRoots\">Polynomial.numberOfRoots</a>.
</p>
</html>"));
  end roots;

  encapsulated function numberOfRoots "Determine number of roots of polynomial"
    import Modelica_LinearSystems2.Math.Polynomial;
    input Polynomial p "Polynomial";
    output Integer result "Number of roots of p";
  protected
    Integer nc=size(p.c, 1);
    Integer i;
  algorithm
    if nc <= 1 then
      // no roots
      result := 0;
    else
      // Remove all leading zero coefficients
      i := 1;
      result := 0;
      while i <= nc - 1 loop
        if p.c[i] <> 0.0 then
          result := nc - i;
          i := nc;
        else
          i := i + 1;
        end if;
      end while;
    end if;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Polynomial.<strong>numberOfRoots</strong>(p);
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate number of roots of a&nbsp;polynomial&nbsp;<code>p</code>.
</p>

<h4>Example</h4>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Polynomial;

x = Polynomial.x();
p1 = -6*x^2 + 4*x -3;

Polynomial.numberOfRoots(p1);
// = 2
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.roots\">Polynomial.roots</a>.
</p>
</html>
"));
  end numberOfRoots;

  encapsulated function rootsOfNonZeroHighestCoefficientPolynomial
    "Determine zeros of polynomial where highest coefficient of polynomial is not zero"
    import Modelica.ComplexMath.j;
    import Modelica_LinearSystems2.Math.Matrices;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Complex;

    input Polynomial p "Polynomial";
    input Integer numberOfRoots "Number of roots of polynomial";
    output Complex result[:]=fill(Complex(0, 0), numberOfRoots)
      "Roots of polynomial";
  protected
    Integer nc=size(p.c, 1);
    Integer i_start=nc - numberOfRoots;
    Integer n=numberOfRoots;
    Real A[n, n] "Companion matrix";
    Real ev[n, 2] "Eigen values";

  algorithm
    assert(numberOfRoots >= 0 and numberOfRoots < nc,
      "Argument numberOfRoots (= " + String(numberOfRoots) +
      ") is not in the range\n" + "0 <= numberOfRoots <= " + String(nc - 1));
    assert(p.c[i_start] <> 0, "p.c[" + String(i_start) +
      "] = 0. Probably wrong argument numberOfRoots (=" + String(numberOfRoots)
       + ")");

    if numberOfRoots > 0 then
      // companion matrix
      A[1, :] := -p.c[i_start + 1:nc]/p.c[i_start];
      A[2:n, :] := [identity(n - 1), zeros(n - 1)];

      // roots are eigenvalues of companion matrix
      //    ev := Matrices.eigenValues(A);
      (ev[:, 1],ev[:, 2]) := Matrices.Internal.eigenvaluesHessenberg(A);
      for i in 1:n loop
        result[i] := ev[i, 1] + j*ev[i, 2];
      end for;
    end if;
  end rootsOfNonZeroHighestCoefficientPolynomial;

  encapsulated package Internal
    "Internal utility functions of record Polynomial that should not be directly utilized by the user"
    extends Modelica.Icons.InternalPackage;
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;

    encapsulated function evaluate_der
      "Evaluate derivative of polynomial at a given abszissa value"
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial p "Polynomial";
      input Real x "Abszissa value";
      input Real dx "Derivative of abszissa value, der(x)";
      output Real dy "Derivative value of polynomial at x";
    protected
      Integer n=size(p.c, 1);

    algorithm
      dy := p.c[1]*(n - 1);
      for j in 2:n - 1 loop
        dy := p.c[j]*(n - j) + x*dy;
      end for;
      dy := dy*dx;

      annotation (
        derivative(order=2)=Polynomial.Internal.evaluate_dder,
        Documentation(info="<html>
<p>
This function is the <em>first time derivative</em> of the function
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluate\">Polynomial.evaluate</a>.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.Internal.evaluate_dder\">Internal.evaluate_dder</a>
</p>
</html>"));
    end evaluate_der;

    encapsulated function evaluate_dder "Evaluate 2nd derivative of polynomial at a given abszissa value"
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial p "Polynomial";
      input Real x "Abszissa value";
      input Real dx "Derivative of abszissa value, der(x)";
      input Real ddx "Second derivative of abszissa value, der(dx)";
      output Real ddy "Second derivative value of polynomial at x";
    protected
      Integer n=size(p.c, 1);
      Real dy "Derivative value of polynomial at x";

    algorithm
      dy := p.c[1]*(n - 1);
      for j in 2:n - 1 loop
        dy := p.c[j]*(n - j) + x*dy;
      end for;

      ddy := p.c[1]*(n - 1)*(n - 2);
      for j in 2:n - 2 loop
        ddy := p.c[j]*(n - j)*(n - j - 1) + x*ddy;
      end for;

      ddy := ddy*dx*dx + dy*ddx;

      annotation (Documentation(info="<html>
<p>
This function is the <em>second time derivative</em> of the function
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluate\">Polynomial.evaluate</a>.
</p>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.Internal.evaluate_der\">Internal.evaluate_der</a>
</p>
</html>"));
    end evaluate_dder;

    encapsulated function integralValue_der
      "Evaluate derivative of integral of polynomial p(x) from x_low to x_high, assuming only x_high as time-dependent (Leibniz rule)"
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial p "Polynomial";
      input Real x_high "High integrand value";
      input Real x_low=0 "Low integrand value, default 0";
      input Real dx_high "High integrand value";
      input Real dx_low=0 "Low integrand value, default 0";
      output Real dintegral=0.0 "Integral of polynomial p from u_low to u_high";
    algorithm
      dintegral := Polynomial.evaluate(p, x_high)*dx_high;
    end integralValue_der;

    function mult
      "Multiply two polynomials (polynomials are defined by vectors)"
      import Modelica.Utilities.Streams.print;

      input Real p1[:];
      input Integer n1
        "Number of coefficients of p1 to be used, i.e., (end-n1+1:end)";
      input Real p2[:];
      input Integer n3_max "Dimension of output vector";
      output Real p3[n3_max];
    protected
      Integer n1_max=size(p1, 1);
      Integer n2=size(p2, 1);
      Integer n3=n1 + n2 - 1;
      Real ck;
    algorithm
      for k in 1:n3 loop
        ck := 0.0;
        for j in max(1, k + 1 - n2):min(k, n1) loop
          ck := ck + p1[n1_max - n1 + j]*p2[k + 1 - j];
        end for;
        p3[n3_max - n3 + k] := ck;
      end for;
    end mult;
  end Internal;

  annotation (defaultComponentName="polynomial", Documentation(info="<html>
<p>
This record defines a polynomial, e.g.,
<var>y</var>&nbsp;=&nbsp;2*<var>x</var><sup>2</sup>&nbsp;+ 3*<var>x</var>&nbsp;+&nbsp;1.
The general form is:
</p>
<pre>
   y = c[1]*x^n + c[2]*x^(n-1) + ... + c[n]*x + c[n+1];
</pre>
<p>
In the record, the coefficients <code>c[i]</code> are stored. Usually,
the record is not directly accessed. Instead, a&nbsp;polynomial is
generated with the functions provided in the record, see
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.'constructor'\">Polynomial.&apos;constructor&apos;</a>.
Also
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.fitting\">Polynomial.fitting(..)</a>
can be helpful.
Several functions are provided that operate on polynomials.
</p>

<h4>
See also
</h4>
<p>
Introduction of the Polynomials in
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.GettingStarted.Polynomials\">Getting started</a> section.
</p>
</html>"));
end Polynomial;
