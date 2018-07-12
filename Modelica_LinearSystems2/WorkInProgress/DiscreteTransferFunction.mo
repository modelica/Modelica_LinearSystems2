within Modelica_LinearSystems2.WorkInProgress;
record DiscreteTransferFunction
  "Discrete transfer function description of a single input, single output system (data + operations)"

  extends Modelica.Icons.Record;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));

  Modelica.SIunits.Time Ts "Sample time"
       annotation(Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

  String uName="u" "Name of input signal"    annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal"  annotation(Dialog(group="Signal names"));

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
      input Modelica.SIunits.Time Ts "Sample time";
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
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Complex;

      input Complex z[:]=fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
      input Complex p[:]=fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
      input Real k=1.0 "Constant multiplied with transfer function";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(
        redeclare Real n[size(z, 1)+1],
        redeclare Real d[size(p, 1)+1])
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
   <b>import</b> Modelica_LinearSystems2.Math.Complex;
   <b>import</b> Modelica_LinearSystems2.ZerosAndPoles;

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
           input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

           input String uName = "" "input name";
           input String yName = "" "output name";

           output DiscreteTransferFunction dtf(redeclare Real n[size(n, 1)], redeclare Real
               d[                                                                             size(d, 1)])
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
      input Modelica.SIunits.Time Ts "Sample time";
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
      input Modelica.SIunits.Time Ts "Sample time"
           annotation(Dialog(group="Data used to construct discrete from continuous system"));

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
z = DiscreteTransferFunction.<b>s</b>()
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

encapsulated package Plot "Functions to plot state space system responses"

encapsulated function bode "Plot transfer function as bode plot"
      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
      import Complex;
      import Modelica_LinearSystems2.Utilities.Plot;
      import SI = Modelica.SIunits;

  input DiscreteTransferFunction dtf "DiscreteTransfer function to be plotted";
  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
        "True, if abszissa range is automatically determined";
  input SI.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false"
                                                    annotation(Dialog(enable=not autoRange));
  input SI.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                              annotation(Dialog(enable=not autoRange));

  input Boolean magnitude=true "= true, to plot the magnitude of tf"
    annotation(choices(checkBox=true));
  input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot of  dtf = "
         + String(dtf)));

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
  TransferFunction tf=TransferFunction(n=dtf.n, d=dtf.d);

  Plot.Records.Curve curves[2];
  Integer i;
  Plot.Records.Diagram diagram2[2];

algorithm
        // Determine frequency vector f
  if autoRange then
    (numZerosZ,denZerosZ) := TransferFunction.Analysis.zerosAndPoles(tf);
  else
    numZerosZ := fill(Complex(0), 0);
    denZerosZ := fill(Complex(0), 0);
  end if;

  numZeros := fill(Complex(0),0);
 // numZeros := fill(Complex(0),size(numZerosZ,1));
  // for i in 1:size(numZerosZ,1) loop
  //   numZeros[i] := Complex.log(numZerosZ[i])/dtf.Ts;
  // end for;

  denZeros := fill(Complex(0),size(denZerosZ,1));
  for i in 1:size(denZerosZ,1) loop
    denZeros[i] :=Modelica.ComplexMath.log(denZerosZ[i])/dtf.Ts;
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
    z[i] :=Modelica.ComplexMath.exp(Complex(0, w[i]*dtf.Ts));
    c := TransferFunction.Analysis.evaluate(
          tf,
          z[i],
          1e-10);
    A[i] :=Modelica.ComplexMath.'abs'(c);
    phi_old :=Modelica.ComplexMath.arg(c, phi_old);
    phi[i] := SI.Conversions.to_deg(phi_old);

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
    diagram2[i].yLabel := "magnitude";
    if phase then
       diagram2[i].xLabel:="";
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

    if magnitude and phase then
      Plot.diagramVector(diagram2, device);
    else
      Plot.diagram(diagram2[1], device);
    end if;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
TransferFunction.Plot.<b>plotBode</b>(dtf)
   or
TransferFunction.Plot.<b>plotBode</b>(dtf, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the bode-diagram of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf =(s^2 + 5*s + 7)/(s^2 + 5*s + 6);

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.plotBode(tf)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</blockquote>
</html>"));
end bode;

end Plot;

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
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Analysis.<b>denominatorDegree</b>(dtf)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of a discrete transfer function.

</p>

<h4>Example</h4>
<blockquote><pre>
   TransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
   Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^2+z+1);

   Real dDegree;

<b>algorithm</b>
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
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteZerosAndPoles.Internal;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteTransferFunction;
  import Modelica_LinearSystems2.TransferFunction;
  import Complex;

  input DiscreteTransferFunction dtf "transfer function of a system";
  output DiscreteZerosAndPoles dzp(
    redeclare Real n1[Internal.numberOfRealZeros2(dtf)],
    redeclare Real n2[integer((size(dtf.n, 1) - 1 -
      Internal.numberOfRealZeros2(dtf))/2),2],
    redeclare Real d1[Internal.numberOfRealPoles(dtf)],
    redeclare Real d2[integer((size(dtf.d, 1) - 1 -
      Internal.numberOfRealPoles(dtf))/2),2]);
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
<table>
<tr> <td align=right>  dzp </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Conversion.<b>toDiscreteZerosAndPoles</b>(tf)  </td> </tr>
</table>

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
of a discrete transfer function representated by numerator and denominator polynomial. The poles and zeros and the gain <tt>k</tt> are computed by
(<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs the DiscreteZerosAndPoles constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.TransferFunction tf = 1/(z^2 + 3*z +2)

<b>algorithm</b>
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
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Conversion.toStateSpace<b>toDiscreteStateSpace</b>(dtf)  </td> </tr>
</table>

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
<b>der</b>(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
   with
           <b>A</b> = [   0  ,    1  ,    0  ,    0;
                   0  ,    0  ,    1  ,    0:
                   0  ,    0  ,    0  ,    1;
                -a0/a4, -a1/a4, -a2/a4, -a3/a4];
            <b>B</b> = [  0;
                  0;
                  0;
                 1/a4];
           <b>C</b> = [b0-b4*a0/a4, b1-b4*a1/a4, b2-b4*a2/a4, b3-b4*a3/a4];
           <b>D</b> = [b4/a4];
</pre></blockquote>
<p>
If the numerator polynomial is 1, then the state vector
<b>x</b> is built up of the y(k) (the privious y) and of all the nx-1 predecessor
(nx is the dimension of the state vector):
</p>
<blockquote><pre>
   <b>x</b>(k+1) = {y(k-n+1), y(k-n+2), ..., y(k)};
</pre></blockquote>
<p>
Note, the state vector <b>x</b> of Modelica.Blocks.Continuous.TransferFunction
is defined slightly differently.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^3 + z^2 + z +1);

<b>algorithm</b>
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
