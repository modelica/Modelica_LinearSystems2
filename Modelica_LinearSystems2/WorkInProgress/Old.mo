within Modelica_LinearSystems2.WorkInProgress;
package Old
  "Old implementations, kept for comparison, if bugs occur in the future"

  encapsulated function ZerosAndPolesEvaluate
    "Evaluate a ZerosAndPoles transfer function at a given value of p"
    import Modelica.Utilities.Streams.print;
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Complex;

    input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
    input Complex p=Complex(0) "Complex value p where zp is evaluated";
    input Real den_min(min=0)=0 "|denominator(p)| is limited by den_min";
    output Complex y "= zp(p)";
  protected
    Complex j=Modelica.ComplexMath.j;
    Complex num;
    Complex den;
    Real abs_den;
  algorithm
    // Build numerator
    num := zp.k + 0*j;
    for i in 1:size(zp.n1, 1) loop
      num := num*ZerosAndPoles.Internal.'p+a'(p, zp.n1[i]);
    end for;
    for i in 1:size(zp.n2, 1) loop
      num := num*ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(p, zp.n2[i, :]);
      print("... evaluate 1: i="+ String(i)+", num = "+ String(num));
    end for;

    // Build denominator
    den := 1 + 0*j;
    for i in 1:size(zp.d1, 1) loop
      den := den*ZerosAndPoles.Internal.'p+a'(p, zp.d1[i]);
    end for;
    for i in 1:size(zp.d2, 1) loop
      den := den*ZerosAndPoles.Internal.'p^2+k[1]*p+k[2]'(p, zp.d2[i, :]);
      print("... evaluate 2: i="+ String(i)+", den = "+ String(den));
    end for;

    // Build value of transfer function
    abs_den :=Modelica.ComplexMath.'abs'(den);
    den := if abs_den >= den_min then den else (if den.re >= 0 then den_min else -
      den_min) + 0*j;
    print("... evaluate 3: num = "+ String(num) + ", den = "+String(den));

    y := num/den;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = ZerosAndPoles.Analysis.<b>evaluate</b>(zp, p, den_min=0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the ZerosAndPoles transfer function at a given (complex) value of p and returns the value G(p)=N(p)/D(p). The optional argument den_min with default 0 is used to guard against a division by zero.
</p>
<pre>
  <b>if</b> |(D(p))| >= den_min <b>then</b>
     G(p) = N(p) / D(p);
  <b>elseif</b> D(p).re >= 0.0 <b>then</b>
     G(p) = N(p) / den_min
  <b>else</b>
     G(p) = -N(p) / den_min
  <b>end if</b>;
</p>
</pre>

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
</html>",   revisions="<html>
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
  end ZerosAndPolesEvaluate;

  function analysisInvariantZeros
    "Example to compute the invariant zeros of a state space system"
    import Modelica;
    import Modelica.Utilities.Streams.print;
    import Complex;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2;

    input Complex z[:]={-2 + 0*j,-3 + 4*j,-3 - 4*j} "Zeros (Complex vector of numerator zeros)";
    input Complex p[:]={-0.5 + 0*j,-5 + 2*j,-5 - 2*j} "Poles (Complex vector of denominator zeros)";
    input Real k=1.0 "Constant multiplied with transfer function";

  protected
    input Complex j = Modelica.ComplexMath.j;

    ZerosAndPoles zp=ZerosAndPoles(
        z=z,
        p=p,
        k=k);

    StateSpace ss=StateSpace(zp);
    Real Poles[:,2];
    Real Zeros[:,2];
    Boolean ok;
  algorithm
    Zeros := StateSpace.Internal.invariantZerosWithRealMatrix(ss.A, ss.B, ss.C, ss.D);
    if size(Zeros, 1) == 0 then
      print("\nSystem\n  "+String(zp)+"\nhas no invariant zeros\n");
    else
      print("\nSystem\n  zp = "+String(zp)+"\n has " + String(size(Zeros, 1)) + " invariant zeros:");
      for i in 1:size(Zeros, 1) loop
        print("   " + String(Zeros[i,1]) + " + j*" + String(Zeros[i,2]));
      end for;
    end if;
    ok := true;

   annotation (Documentation(info="<html>
<p>
This example shows the computation of the poles and zeros of state space system.
</p>
</html>"));
  end analysisInvariantZeros;

  function analysisZerosAndPoles
    "Example to compute the invariant zeros of a state space system"
    import Modelica;
    import Modelica.Utilities.Streams.print;
    import Complex;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2;

    input Complex z[:]={-2 + 0*j,-3 + 4*j,-3 - 4*j} "Zeros (Complex vector of numerator zeros)";
    input Complex p[:]={-0.5 + 0*j,-5 + 2*j,-5 - 2*j} "Poles (Complex vector of denominator zeros)";
    input Real k=1.0 "Constant multiplied with transfer function";

  protected
    input Complex j = Modelica.ComplexMath.j;

    ZerosAndPoles zp=ZerosAndPoles(
        z=z,
        p=p,
        k=k);

    StateSpace ss=StateSpace(zp);
    Real Poles[:,2];
    Real Zeros[:,2];
    Boolean ok;
  algorithm
    (Poles,Zeros) :=
      Modelica_LinearSystems2.WorkInProgress.polesAndZerosAsRealMatrix(
        ss.A,
        ss.B,
        ss.C,
        ss.D);
    print("Poles = " + Modelica.Math.Matrices.toString(Poles));
    print("Zeros = " + Modelica.Math.Matrices.toString(Zeros));

   annotation (Documentation(info="<html>
<p>
This example shows the computation of the poles and zeros of state space system.
</p>
</html>"));
  end analysisZerosAndPoles;
end Old;
