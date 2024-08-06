within Modelica_LinearSystems2.WorkInProgress;
record DiscreteTransferFunction
  "Discrete transfer function description of a single input, single output system (data + operations)"

  extends Modelica.Icons.Record;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));

  Modelica.Units.SI.Time Ts "Sample time"
    annotation (Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

  String uName="u" "Name of input signal" annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal" annotation(Dialog(group="Signal names"));

/* If the numerator polynomial has no coefficients, the transfer function
   is zero. The denominator polynomial must always have at
   least one coefficient, such as {1}
*/

  encapsulated operator 'constructor'
    "Default constructor for a transfer function"
    import Modelica;
    import Modelica_LinearSystems2.TransferFunction;

    function fromReal
      "Generate a DiscreteTransferFunction data record from a Real value"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;

      input Real r "Value of Real variable";
      input Modelica.Units.SI.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n={r}, d={1});

    algorithm
      dtf.Ts := Ts;
      dtf.method := method;
      dtf.uName := uName;
      dtf.yName := yName;
    end fromReal;

    encapsulated function fromZerosAndPoles
      "Generate a discrete transfer function from a set of zeros and poles"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Complex z[:]=fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
      input Complex p[:]=fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
      input Real k=1.0 "Constant multiplied with transfer function";
      input Modelica.Units.SI.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(redeclare Real n[size(z, 1)+1], redeclare Real
               d[size(p, 1)+1])
        "TransferFunction built by ZerosAndPoles object";

    protected
      Polynomial pn=k*Polynomial(z);
      Polynomial pd=Polynomial(p);
    algorithm

      dtf.n := pn.c;
      dtf.d := pd.c;
      dtf.Ts := Ts;
      dtf.method := method;
      dtf.uName := uName;
      dtf.yName := yName;

      annotation (Documentation(info="<html>
<p>
This function constructs a transfer function from numerator zeros.
Example:
</p>

<p> The transfer function</p>
<pre>
   zp = (s+(2+3*j))*(s+(2-3*j))
</pre>
<p>
can be expressed as
</p>
<pre>
   <strong>import</strong> Complex;
   <strong>import</strong> Modelica_LinearSystems2.ZerosAndPoles;

   j = Complex.j();
   zp = ZerosAndPoles({2+3*j}, {2-3*j});
</pre>

<p>
Since only transfer functions with real coefficients are supported,
complex zeros must be defined as conjugate complex pairs.
It is required that complex conjugate pairs must directly
follow each other as above. An error occurs if this is not the case.
</p>
</html>"));
    end fromZerosAndPoles;

    encapsulated function fromArrays
      "Generate a DiscreteTransferFunction data record from numerator and denominator array"
      import Modelica;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2;

      input Real n[:] "Coefficients of numerator polynomial";
      input Real d[:] "Coefficients of denominator polynomial";
      input Modelica.Units.SI.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

      input String uName = "" "input name";
      input String yName = "" "output name";

      output DiscreteTransferFunction dtf(redeclare Real n[size(n, 1)], redeclare Real d[size(d, 1)])
        "Transfer function";

    algorithm
      //this is the constructor algorithm
      assert(size(d, 1) > 0, "Input denominator d must have at least one element, however\n"
         + "d is an empty vector");
      dtf.n := n;
      dtf.d := d;
      dtf.Ts := Ts;
      dtf.method := method;
      dtf.uName := uName;
      dtf.yName := yName;
    end fromArrays;

    function fromPolynomials
      "Generate a DiscreteTransferFunction data record from a numerator and denominator polynomial"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial n "Numerator polynomial";
      input Polynomial d "Denominator polynomial";
      input Modelica.Units.SI.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n=n.c, d=d.c, Ts=Ts, method=method, uName=uName, yName=yName);

    algorithm
    end fromPolynomials;

    function fromTransferFunction
      "Generate a DiscreteTransferFunction data record from a continuous Transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

      input TransferFunction tf "continuous transfer function";
      input Modelica.Units.SI.Time Ts "Sample time" annotation (Dialog(group=
              "Data used to construct discrete from continuous system"));

      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

      output DiscreteTransferFunction dtf;
    protected
      StateSpace ss = StateSpace(tf);
      Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace dss=
                               Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace(
                                                  ss,Ts,method);

    algorithm
      dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss);
    end fromTransferFunction;
  end 'constructor';

    encapsulated operator function 'String'
    "Transform TransferFunction into a String representation"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;

    input DiscreteTransferFunction dtf
      "Discrete transfer function to be transformed in a String representation";
    input Integer significantDigits=6
      "Number of significant digits that are shown";
    input String name="z" "Independent variable name used for printing";
    output String z="";
  protected
    Integer n_num=size(dtf.n, 1) - 1;
    Integer n_den=size(dtf.d, 1) - 1;
    Boolean numParenthesis;
    algorithm
    if n_num == -1 then
      z := "0";
    else
      numParenthesis := n_num > 0 and not (n_den == 0 and dtf.d[1] == 1);
      if numParenthesis then
        z := "(";
      end if;
       z := z + String(
            Polynomial(dtf.n),
            significantDigits,
            name);

      if numParenthesis then
        z := z + ")";
      end if;
    if n_den > 0 or dtf.d[1] <> 1 then
        if n_den > 0 then
          z := z + "/(";
        else
          z := z + "/";
        end if;

        z := z + String(
              Polynomial(dtf.d),
              significantDigits,
              name);

        if n_den > 0 then
          z := z + ")";
        end if;
      end if;
    end if;
    end 'String';

encapsulated function z "Generate the discrete transfer function z"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

  output DiscreteTransferFunction dtf(n={1,0}, d={1}) "z";
algorithm

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
z = DiscreteTransferFunction.<strong>s</strong>()
</pre></blockquote>

<h4>Description</h4>
<p>
Generate the complex variable z=exp(T*s) as a DiscreteTransferFunction. It can be used for generating like
</p>
<blockquote><pre>
DiscreteTransferFunction dtf = z/(3*z^2 + 2*z +2)
</pre></blockquote>
</html>"));
end z;

encapsulated package Analysis

encapsulated function denominatorDegree "Return denominator degree"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;

  input DiscreteTransferFunction dtf "discrete transfer function of a system";
  output Integer result;

algorithm
  result := size(dtf.d,1)-1;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = DiscreteTransferFunction.Analysis.<strong>denominatorDegree</strong>(dtf)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<strong>denominatorDegree</strong> calculates the degree of the denominator polynomial of a discrete transfer function.

</p>

<h4>Example</h4>
<blockquote><pre>
   TransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
   Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^2+z+1);

   Real dDegree;

<strong>algorithm</strong>
  dDegree := TransferFunction.Analysis.denominatorDegree(dtf);
//  dDegree = 2
</pre></blockquote>


</html>"));
end denominatorDegree;

end Analysis;

encapsulated package Conversion

encapsulated function toDiscreteZerosAndPoles
  "Generate a DiscreteZerosAndPoles object from a DiscreteTransferFunction object"
  import Modelica;
  import Complex;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
  import Modelica_LinearSystems2.TransferFunction;

  input DiscreteTransferFunction dtf "transfer function of a system";
  output DiscreteZerosAndPoles dzp(
    redeclare Real n1[DiscreteZerosAndPoles.Internal.numberOfRealZeros2(dtf)],
    redeclare Real n2[integer((size(dtf.n, 1) - 1 -
      DiscreteZerosAndPoles.Internal.numberOfRealZeros2(dtf))/2),2],
    redeclare Real d1[DiscreteZerosAndPoles.Internal.numberOfRealPoles(dtf)],
    redeclare Real d2[integer((size(dtf.d, 1) - 1 -
      DiscreteZerosAndPoles.Internal.numberOfRealPoles(dtf))/2),2]);
    protected
  TransferFunction tf=TransferFunction(n=dtf.n, d=dtf.d);
  Complex z[:];
  Complex p[:];
  Real k;
algorithm
  (z,p,k) := TransferFunction.Analysis.zerosAndPoles(tf);
  dzp := DiscreteZerosAndPoles(z, p, k, dtf.Ts, dtf.method, uName=dtf.uName, yName=dtf.yName);
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteTransferFunction.Conversion.<strong>toDiscreteZerosAndPoles</strong>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a DiscreteZerosAndPoles record
</p>
<blockquote><pre>
          product(z + n1[i]) * product(z^2 + n2[i,1]*z + n2[i,2])
zp = k * ---------------------------------------------------------
          product(z + d1[i]) * product(z^2 + d2[i,1]*z + d2[i,2])
</pre></blockquote>
<p>
of a discrete transfer function represented by numerator and denominator polynomial. The poles and zeros and the gain <tt>k</tt> are computed by
(<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs the DiscreteZerosAndPoles constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.TransferFunction tf = 1/(z^2 + 3*z +2)

<strong>algorithm</strong>
  dzp:=Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles(dtf);
//  zp = 1/( (z + 1)*(z + 2) )
</pre></blockquote>
</html>"));
end toDiscreteZerosAndPoles;

function toDiscreteStateSpace
      "Convert a DiscreteTransferFunction into a DiscreteStateSpace representation"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;
      import Modelica.Math.Vectors;

 input DiscreteTransferFunction dtf "discrete transfer function of a system";
      output DiscreteStateSpace dss(
        redeclare Real A[DiscreteTransferFunction.Analysis.denominatorDegree(dtf),DiscreteTransferFunction.Analysis.denominatorDegree(dtf)],
        redeclare Real B[DiscreteTransferFunction.Analysis.denominatorDegree(dtf),1],
        redeclare Real B2[DiscreteTransferFunction.Analysis.denominatorDegree(dtf),1],
        redeclare Real C[1,DiscreteTransferFunction.Analysis.denominatorDegree(dtf)],
        redeclare Real D[1,1]) "Discrete state space record";

    protected
 Integer na=DiscreteTransferFunction.Analysis.denominatorDegree(dtf) + 1;
 Integer nb=size(dtf.n,1);//numerator degree
 Integer nx=na - 1;
 TransferFunction tf=TransferFunction(n=dtf.n, d=dtf.d);
 Real a[na]=Vectors.reverse(tf.d) "Reverse element order of tf.a";
 Real b[na]=vector([Vectors.reverse(tf.n); zeros(na - nb, 1)]);
 Real d=b[na]/a[na];
algorithm
 if nx == 0 then
   dss.A := fill(0, 0, nx);
   dss.B := fill(0, 0, 1);
   dss.B2 := fill(0, 0, 1);
   dss.C := fill(0, 1, 0);
 else
   dss.A[1:nx - 1, 1:nx] := [zeros(nx - 1, 1),identity(nx - 1)];
   dss.A[nx, 1:nx] := -a[1:na - 1]/a[na];
   dss.B := [zeros(nx - 1, 1); 1/a[na]];
   dss.B2 := fill(0,nx,1);
   dss.C := {b[1:nx]};

end if;
  dss.D := [d];
  dss.Ts := dtf.Ts;
  dss.method := dtf.method;

 annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = DiscreteTransferFunction.Conversion.toStateSpace<strong>toDiscreteStateSpace</strong>(dtf)
</pre></blockquote>

<h4>Description</h4>
<p>
Transforms a discrete transfer function into discrete state space representation.
There are an infinite number of possible realizations.
Here, the transfer function is transformed into
controller canonical form, i.e. the transfer function
</p>
<blockquote><pre>
     b4*z^4 + b3*z^3 + b2*z^2 + b1*z + b0
y = -------------------------------------- *u
     a4*z^4 + a3*z^3 + a2*z^2 + a1*z + a0
</pre></blockquote>
<p>
is transformed into:
</p>
<blockquote><pre>
<strong>der</strong>(<strong>x</strong>) = <strong>A</strong>*<strong>x</strong> + <strong>B</strong>*<strong>u</strong>;
    <strong>y</strong>  = <strong>C</strong>*<strong>x</strong> + <strong>D</strong>*<strong>u</strong>;
   with
           <strong>A</strong> = [   0  ,    1  ,    0  ,    0;
                   0  ,    0  ,    1  ,    0:
                   0  ,    0  ,    0  ,    1;
                -a0/a4, -a1/a4, -a2/a4, -a3/a4];
            <strong>B</strong> = [  0;
                  0;
                  0;
                 1/a4];
           <strong>C</strong> = [b0-b4*a0/a4, b1-b4*a1/a4, b2-b4*a2/a4, b3-b4*a3/a4];
           <strong>D</strong> = [b4/a4];
</pre></blockquote>
<p>
If the numerator polynomial is 1, then the state vector
<strong>x</strong> is built up of the y(k) (the previous y) and of all the nx-1 predecessor
(nx is the dimension of the state vector):
</p>
<blockquote><pre>
   <strong>x</strong>(k+1) = {y(k-n+1), y(k-n+2), ..., y(k)};
</pre></blockquote>
<p>
Note, the state vector <strong>x</strong> of Modelica.Blocks.Continuous.TransferFunction
is defined slightly differently.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^3 + z^2 + z +1);

<strong>algorithm</strong>
  dss := Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteStateSpace(dtf);
// dss.A = [0, 1, 0; 0, 0, 1; -1, -1, -1],
// dss.B = [0; 0; 1],
// dss.C = [1, 1, 0],
// dss.D = [0],
</pre></blockquote>

</html>"));
end toDiscreteStateSpace;

end Conversion;

end DiscreteTransferFunction;
