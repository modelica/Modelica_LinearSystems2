within Modelica_LinearSystems2;
operator record DiscreteTransferFunction
  "Discrete transfer function description of a single input, single output system (data + operations)"

  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{z^m, ... , z, 1} / (d*{z^r, ... , z, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{z^m, ... , z, 1} / (d*{z^r, ... , z, 1}) * u"));

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
    "Collection of operators to construct a DiscreteTransferFunction data record"
    import Modelica;
    import Modelica_LinearSystems2.TransferFunction;

    encapsulated function fromReal
      "Generate a DiscreteTransferFunction data record from a real value"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts=1 "Sample time";
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
      "Generate a DiscreteStateSpace data record from a set of zeros and poles"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Complex z[:] = fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
      input Complex p[:] = fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
      input Real k=1.0 "Constant multiplied with transfer function";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(redeclare Real n[size(z, 1)+1], redeclare Real
               d[                                                                          size(p, 1)+1])
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
</html>"));
    end fromZerosAndPoles;

    encapsulated function fromArrays
      "Generate a DiscreteTransferFunction data record from numerator and denominator array"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2;

      input Real n[:] "Coefficients of numerator polynomial";
      input Real d[:] "Coefficients of denominator polynomial";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

      input String uName = "" "Input name";
      input String yName = "" "Output name";
      output DiscreteTransferFunction dtf(redeclare Real n[size(n, 1)], redeclare Real d[size(d, 1)])
        "Transfer function";

    algorithm
      // This is the constructor algorithm
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
      import Modelica_LinearSystems2.DiscreteTransferFunction;
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
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input TransferFunction tf "continuous transfer function";
      input Modelica.SIunits.Time Ts "Sample time"
           annotation(Dialog(group="Data used to construct discrete from continuous system"));

      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

      output DiscreteTransferFunction dtf;
    protected
      StateSpace ss = StateSpace(tf);
      DiscreteStateSpace dss=DiscreteStateSpace(ss,Ts,method);

    algorithm
      dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss);
    end fromTransferFunction;
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

  encapsulated operator '+'
    "Contains operators for addition of discrete transfer functions"
    import Modelica;

    function 'dtf+dtf'
      "Parallel connection of two discrete transfer functions (= inputs are the same, outputs of the two systems are added)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf1 "Transfer function system 1";
      input DiscreteTransferFunction dtf2 "Transfer function system 2";
      output DiscreteTransferFunction result;

    algorithm
      assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps,
        "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"+\".");
      result := DiscreteTransferFunction(Polynomial(dtf1.n)*Polynomial(dtf2.d) + Polynomial(dtf2.n)*Polynomial(dtf1.d), Polynomial(dtf1.d)*Polynomial(dtf2.d), Ts=dtf1.Ts, method=dtf1.method);
    end 'dtf+dtf';

    function 'dtf+r' "Add a real number to a DiscreteTransferFunction"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf "Transfer function system";
      input Real r "Real number";
      output DiscreteTransferFunction result;
    protected
      DiscreteTransferFunction dtfr=DiscreteTransferFunction(r, Ts=dtf.Ts, method=dtf.method);
    algorithm
      result := dtf + dtfr;
    end 'dtf+r';
    annotation (Documentation(info="<html>
<p>This package contains operators for addition of discrete fransfer function records. </p>
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
            points={{0,50},{0,-50}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end '+';

  encapsulated operator '-'
    "Contains operators for subtraction of discrete transfer functions"
    import Modelica;

    function subtract "Subtract two discrete transfer functions (dtf1 - dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf1;
      input DiscreteTransferFunction dtf2;

      output DiscreteTransferFunction result;

    protected
      Polynomial n=Polynomial(dtf1.n)*Polynomial(dtf2.d) - Polynomial(dtf2.n)*Polynomial(dtf1.d);
      Polynomial d=Polynomial(dtf1.d)*Polynomial(dtf2.d);

    algorithm
      if size(dtf1.n,1)>1 and size(dtf2.n,1)>1 then
        assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"-.subtract\".");
      end if;
      result := DiscreteTransferFunction(n=n.c, d=d.c, Ts=dtf1.Ts, method=dtf1.method);
    end subtract;

    function negate "Unary minus (multiply discrete transfer function by -1)"
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf;
      output DiscreteTransferFunction result(n=-dtf.n, d=dtf.d, Ts=dtf.Ts, method=dtf.method) "= -dtf";
    algorithm
    end negate;
    annotation (Documentation(info="<html>
<p>This package contains operators for subtraction of discrete fransfer function records. </p>
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
    "Contains operators for multiplication of discrete transfer functions"
    import Modelica;

    function 'dtf*dtf' "Multiply two discrete transfer functions (dtf1 * dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf1 "Transfer function system 1";
      input DiscreteTransferFunction dtf2 "Transfer function system 2";
      output DiscreteTransferFunction result;

    algorithm
      assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"+\".");
      result := DiscreteTransferFunction(Polynomial(dtf1.n)*Polynomial(dtf2.n),Polynomial(dtf1.d)*Polynomial(dtf2.d), Ts=dtf1.Ts, method=dtf1.method);
    end 'dtf*dtf';

    function 'r*dtf'
      "Multiply a real number with a DiscreteTransferFunctions (r * dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input Real r "Real number";
      input DiscreteTransferFunction dtf "Transfer function system";
      output DiscreteTransferFunction result=dtf;

    algorithm
      result.n := r * dtf.n;
    end 'r*dtf';


    annotation (Documentation(info="<html>
<p>This package contains operators for multiplication of discrete fransfer function records. </p>
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

  encapsulated operator '/'
    "Contains operators for division of discrete transfer functions"
    import Modelica;

    encapsulated function 'dtf/dtf'
      "Divide two discrete transfer functions (dtf1 / dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf1 "Transfer function system 1";
      input DiscreteTransferFunction dtf2 "Transfer function system 2";
      output DiscreteTransferFunction result;

    algorithm
      assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"/\".");
      result := DiscreteTransferFunction(Polynomial(dtf1.n)*Polynomial(dtf2.d),Polynomial(dtf1.d)*
        Polynomial(dtf2.n),Ts=dtf1.Ts, method=dtf1.method);
    end 'dtf/dtf';

    function 'r/dtf'
      "Divide a real number by  discrete transfer functions (r / dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input Real r "Real number";
      input DiscreteTransferFunction dtf "Transfer function system";
      output DiscreteTransferFunction result;

    algorithm
      result := DiscreteTransferFunction(r*dtf.d,dtf.n,Ts=dtf.Ts, method=dtf.method);
    end 'r/dtf';


    annotation (Documentation(info="<html>
<p>This package contains operators for division of discrete fransfer function records. </p>
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
            radius=25.0),
          Line(
            points={{20,50},{-20,-50}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end '/';

  encapsulated operator function '^'
    "Integer power of DiscreteTransferFunction (dtf1^k)"
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.DiscreteTransferFunction;
     import Modelica_LinearSystems2.TransferFunction;

     input DiscreteTransferFunction dtf "Transfer function";
     input Integer k(min=0) = 1 "Integer exponent";
     output DiscreteTransferFunction result;

  protected
    TransferFunction tf=(Polynomial(dtf.n)^k)/(Polynomial(dtf.d)^k);

  algorithm
       result := DiscreteTransferFunction(n=tf.n, d=tf.d, Ts=dtf.Ts, method=dtf.method);
  end '^';

  encapsulated operator function '=='
    "Check whether two discrete transfer functions are identical"

     import Modelica;
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.DiscreteTransferFunction;

     input DiscreteTransferFunction dtf1 "Transfer function system 1";
     input DiscreteTransferFunction dtf2 "Transfer function system 1";
     input Real eps(min=0) = 0
      "Two coefficients c1 and c2 of the two transfer functions are identical if abs(c1-c2) <= eps";

     output Boolean result "= dtf1 == dtf2";
  algorithm
    result := (Polynomial(dtf1.n) == Polynomial(dtf2.n)) and (Polynomial(dtf1.d) == Polynomial(dtf2.d) and abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps);
  end '==';

    encapsulated operator function 'String'
    "Transform DiscreteTransferFunction into a String representation"
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

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
      z := z + "\n\n Ts = " + String(dtf.Ts) + "\n method =" +
        Modelica_LinearSystems2.Internal.methodString(dtf.method);

    end 'String';

  encapsulated function z "Generate the discrete transfer function z"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    input Modelica.SIunits.Time Ts=0;
    output DiscreteTransferFunction dtf(n={1,0}, d={1},Ts=Ts) "z";
  algorithm

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
z = DiscreteTransferFunction.<b>z</b>()
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
    "Package of functions to analyse discrete transfer function represented by a DiscreteTransferFunction record"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function timeResponse
      "Calculate the time response of a discrete transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;// Input/Output declarations of discrete time response functions
      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step;
      input Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

    protected
      DiscreteStateSpace dss=DiscreteStateSpace(dtf);
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
(y, t, x) = DiscreteTransferFunction.Analysis.<b>timeResponse</b>(dtf, tSpan, responseType, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the discrete transfer function representation is transformed into a discrete state space representation which is given to DiscreteStateSpace.Analysis.timeResponse and the time response of the state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step    &quot;Step response&quot;,
Ramp    &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;
</pre></blockquote>
<p>
The outputs y and x are calculated from the system equations of the discrete state space system for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[1]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  dtf.Ts:=0.1;
  (y,t,x):=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(dtf,tSpan,response,x0);
//  y[:,1,1] = {0.00237529691211404, 0.0116282350020595, 0.0293927396867651, 0.0546913271597482, 0.0865678034508828}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0.0, 0.0, 1.0, 2.89548693586698, 5.58336953639396}
</pre></blockquote>
</html>"));
    end timeResponse;

    encapsulated function impulseResponse
      "Calculate the impulse time response of a discrete transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

    protected
      Real tSpanVar;
    algorithm
      // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dtf));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=DiscreteTransferFunction.Analysis.timeResponse(
            dtf=dtf,
            tSpan=tSpanVar,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
            x0=zeros(DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = DiscreteTransferFunction.Analysis.<b>impulseResponse</b>(dtf, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the discrete transfer function representation is transformed into discrete state space representation which is given to DiscreteStateSpace.Analysis.timeResponse
and the impulse response of the discrete state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. in this case
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
</pre></blockquote>
<p>
The outputs y and x of the discrete state space systrem are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real tSpan= 0.4;
  Real x0[1]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  dtf.Ts:=0.1;
  (y,t,x):=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.impulseResponse(dtf,tSpan,response,x0);
//  y[:,1,1] = {0.00237529691211404, 0.00925293808994548, 0.0177645046847056, 0.0252985874729831, 0.0318764762911345}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0.0, 0.0, 1.0, 1.89548693586698, 2.68788260052697}
</pre></blockquote>
</html>"));
    end impulseResponse;

    encapsulated function stepResponse
      "Calculate the step time response of a discrete transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;
    protected
      Real tSpanVar;
    algorithm
      // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dtf));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=DiscreteTransferFunction.Analysis.timeResponse(
            dtf=dtf,
            tSpan=tSpanVar,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
            x0=zeros(DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = DiscreteTransferFunction.Analysis.<b>stepResponse</b>(dtf, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the discrete transfer function representation is transformed into discrete state space representation which is given to DiscreteStateSpace.Analysis.timeResponse
and the step response of the discrete state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. in this case
</p>
<blockquote><pre>
Step &quot;Step response&quot;,
</pre></blockquote>
<p>
The outputs y and x of the discrete state space systrem are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real tSpan= 0.4;
  Real x0[1]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  dtf.Ts:=0.1;
  (y,t,x):=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.stepResponse(dtf,tSpan,response,x0);
//  y[:,1,1] = {0.00237529691211404, 0.0116282350020595, 0.0293927396867651, 0.0546913271597482, 0.0865678034508828}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0.0, 0.0, 1.0, 2.89548693586698, 5.58336953639396}
</pre></blockquote>
</html>"));
    end stepResponse;

    encapsulated function rampResponse
      "Calculate the ramp time response of a discrete transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

    protected
      Real tSpanVar;
    algorithm
      // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dtf));
      else
        tSpanVar := tSpan;
      end if;
      (y,t,x_discrete) :=DiscreteTransferFunction.Analysis.timeResponse(
            dtf=dtf,
            tSpan=tSpanVar,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
            x0=zeros(DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = DiscreteTransferFunction.Analysis.<b>rampResponse</b>(dtf, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the discrete transfer function representation is transformed into discrete state space representation which is given to DiscreteStateSpace.Analysis.timeResponse
and the ramp response of the discrete state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. in this case
</p>
<blockquote><pre>
Ramp &quot;Ramp response&quot;,
</pre></blockquote>
<p>
The outputs y and x of the discrete state space systrem are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real tSpan= 0.4;
  Real x0[1]={0,0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  dtf.Ts:=0.1;
  (y,t,x):=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.rampResponse(dtf,tSpan,response,x0);
//  y[:,1,1] = {0.0, 0.000237529691211404, 0.00140035319141736, 0.00433962716009387, 0.00980875987606869}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0.0, 0.0, 0.0, 0.1, 0.389548693586699}
</pre></blockquote>
</html>"));
    end rampResponse;

    encapsulated function initialResponse
      "Calculate the initial time response of a discrete transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input Real x0[:]=fill(0,0) "Initial state vector";

    // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

    protected
      Real tSpanVar;
    algorithm
      // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dtf));
      else
        tSpanVar := tSpan;
      end if;
      (y,t,x_discrete) :=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(
            dtf=dtf,
            tSpan=tSpanVar,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
            x0=x0);

      annotation (
        Documentation(info="<html>
 <h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = DiscreteTransferFunction.Analysis.<b>initialResponse</b>(x0, dtf, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the discrete transfer function representation is transformed into discrete state space representation which is given to DiscreteStateSpace.Analysis.timeResponse
and the initial response of the discrete state space system for initial state x0 is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. in this case
</p>
<blockquote><pre>
Initial &quot;Initial response&quot;,
</pre></blockquote>
<p>
The outputs y and x of the discrete state space systrem are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real tSpan= 0.4;
  Real x0[1]={1,2};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  dtf.Ts:=0.1;
  (y,t,x):=Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.initialResponse(x0,dtf,tSpan,response,x0);
//  y[:,1,1] = {0.0187315575967189, 0.0271552102903869, 0.0345205091861731, 0.0408580313775029, 0.0462052138701078}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {1.0, 2.0, 2.88598574821853, 3.66037203581564, 4.3264045475288}
</pre></blockquote>
</html>"));
    end initialResponse;

    encapsulated function denominatorDegree
      "Return denominator degree of a discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf
        "discrete transfer function of a system";
      output Integer result;

    algorithm
      result := size(dtf.d,1)-1;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = DiscreteTransferFunction.Analysis.<b>denominatorDegree</b>(dtf)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of a discrete transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction  dtf=(0.0023753*z^2 + 0.00475059*z + 0.0023753)/(z^2 - 1.89549*z + 0.904988);

  Real dDegree;

<b>algorithm</b>
  dDegree := DiscreteTransferFunction.Analysis.denominatorDegree(dtf);
//  dDegree = 2
</pre></blockquote>
</html>"));
    end denominatorDegree;

  end Analysis;

  encapsulated package Plot
    "Package of functions to plot discrete transfer function responses"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function bode "Plot discrete transfer function as bode plot"
      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Plot;
      import Complex;
      import SI = Modelica.SIunits;

      input DiscreteTransferFunction dtf
        "DiscreteTransfer function to be plotted";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange=true
        "True, if abszissa range is automatically determined";
      input SI.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false" annotation(Dialog(enable=not autoRange));
      input SI.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false" annotation(Dialog(enable=not autoRange));

      input Boolean magnitude=true "= true, to plot the magnitude of tf"
                                                                        annotation(choices(checkBox=true));
      input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot: "
             + String(dtf)));

      input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));

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

      annotation (Documentation(info="<html>
</html>"));
    end bode;

    encapsulated function timeResponse
      "Plot the time response of a system represented by a discrete transfer function. The response type is selectable"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input Modelica_LinearSystems2.DiscreteTransferFunction dtf;
    //  input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "Type of time response";
      input Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
          DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Time response of  dtf = "
             + String(dtf)));

    protected
      Plot.Records.Curve curve;
      Plot.Records.Diagram diagram2;
      Real y[:,1,1] "Output response";
      Real t[:] "Time vector: (number of samples)";

      Real yy[:] "Output response";
      Real tt[:] "Time vector: (number of samples)";

    algorithm
      (y,t) := DiscreteTransferFunction.Analysis.timeResponse(
        dtf,
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

      annotation (
        Documentation(info="<html>
</html>"));
    end timeResponse;

    encapsulated function impulse
      "Impulse response plot of a discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Plot;

      input DiscreteTransferFunction dtf "zeros-and-poles transfer function";
      input Real tSpan=0 "Simulation time span [s]";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
           Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response of  zp = "
             + String(dtf)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "Type of time response";

      Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
           DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";
    algorithm
      // set sample time
      Modelica_LinearSystems2.DiscreteTransferFunction.Plot.timeResponse(
          dtf=dtf,
          tSpan=tSpan,
          response=response,
          x0=x0,
          defaultDiagram=defaultDiagram,
          device=device);

      annotation (
        Documentation(info="<html>
</html>"));
    end impulse;

    encapsulated function step
      "Step response plot of a discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

      input DiscreteTransferFunction dtf;
      input Real tSpan=0 "Simulation time span [s]";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  dtf = "
             + String(dtf)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
      Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
          DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

    algorithm
      DiscreteTransferFunction.Plot.timeResponse(
        dtf=dtf,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (
        Documentation(info="<html>
</html>"));
    end step;

    encapsulated function ramp
      "Ramp response plot of a discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input DiscreteTransferFunction dtf;
      input Real tSpan=0 "Simulation time span [s]";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  dtf = "
             + String(dtf)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";
      Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
          DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";
    algorithm
     Modelica_LinearSystems2.DiscreteTransferFunction.Plot.timeResponse(
        dtf=dtf,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (
        Documentation(info="<html>
</html>"));
    end ramp;

    encapsulated function initialResponse
      "Initial condition response plot of a discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input Modelica_LinearSystems2.DiscreteTransferFunction dtf;
      input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";
      input Real y0 "Initial output (for initial condition plot)";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Initial response of  dtf = "
             + String(dtf) + "  with y0 = " + String(y0)));

    protected
      Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(dtf);
      Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=
          Modelica.Math.Matrices.equalityLeastSquares(
          dss.A,
          fill(0, size(dss.B, 1)),
          dss.C,
          vector(y0)) "Initial state vector (for initial condition plot)";
    algorithm
      Modelica_LinearSystems2.DiscreteTransferFunction.Plot.timeResponse(
            dtf=dtf,
            tSpan=tSpan,
            response=response,
            x0=x0,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
</html>"));
    end initialResponse;

  end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of DiscreteTransferFunction data record"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function toDiscreteZerosAndPoles
      "Generate a DiscreteZerosAndPoles object from a DiscreteTransferFunction object"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;

      input DiscreteTransferFunction dtf "transfer function of a system";
      output Modelica_LinearSystems2.DiscreteZerosAndPoles dzp(
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
      dzp := DiscreteZerosAndPoles(
            z,
            p,
            k,
            dtf.Ts,
            dtf.method,
            uName=dtf.uName,
            yName=dtf.yName);
      annotation (Documentation(info=
                                 "<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteTransferFunction.Conversion.<b>toDiscreteZerosAndPoles</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a DiscreteZerosAndPoles record
<blockquote><pre>
           product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
dzp = k * ---------------------------------------------------------
           product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre></blockquote>
<p>
of a discrete transfer function representated by numerator and denominator
polynomial. The poles and zeros and the gain <tt>k</tt> are computed from
the DiscreteTransferFunction-input and are used as inputs the
DiscreteZerosAndPoles constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf = 1/(z^2 + 3*z + 2)


<b>algorithm</b>
  dzp = Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteZerosAndPoles(dtf);
//  dzp = 1/( (z + 1)*(z + 2) )
</pre></blockquote>
</html>"));
    end toDiscreteZerosAndPoles;

    function toDiscreteStateSpace
      "Convert a DiscreteTransferFunction into a DiscreteStateSpace representation"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica.Math.Vectors;

     input DiscreteTransferFunction dtf
        "discrete transfer function of a system";
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
     Real b[na];//=vector([Vectors.reverse(tf.n); zeros(na - nb, 1)]);
     Real d;//=b[na]/a[na];
    algorithm
     assert(nb<=na,"DiscreteTransferFunction\n" +String(dtf) +"\nis acausal and cannot be transformed to DiscreteStaeSpace in function \"Conversion.toDiscreteStateSpace()\"");
     b := vector([Vectors.reverse(tf.n); zeros(na - nb, 1)]);
     d := b[na]/a[na];
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
       dss.C := {b[1:nx] - d*a[1:nx]};

    end if;
      dss.D := [d];
      dss.Ts := dtf.Ts;
      dss.method := dtf.method;

     annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = DiscreteTransferFunction.Conversion<b>toDiscreteStateSpace</b>(dtf)
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
y = -------------------------------------- * u
     a4*z^4 + a3*z^3 + a2*z^2 + a1*z + a0
</pre></blockquote>
<p>
is transformed into:
</p>
<blockquote><pre>
<b>x</b>_k+1 = <b>A</b>*<b>x</b>_k + <b>B</b>*<b>u</b>_k;
 <b>y</b>_k  = <b>C</b>*<b>x</b>_k + <b>D</b>*<b>u</b>_k;
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

<h4>Example</h4>
<blockquote><pre>
  DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
  Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^3 + z^2 + z + 1);

<b>algorithm</b>
  dss := Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteStateSpace(dtf);
// ss.A = [0, 1, 0; 0, 0, 1; -1, -1, -1],
// ss.B = [0; 0; 1],
// ss.C = [1, 1, 0],
// ss.D = [0],
// ss.B2 = [0; 0; 0],
</pre></blockquote>
</html>"));
    end toDiscreteStateSpace;

  end Conversion;

  encapsulated package Import
    "Package of functions to generate a DiscreteTransferFunction data record from imported data"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function fromFile
      "Generate a DiscreteTransferFunction data record by reading numenator coefficients and denominator coefficients from a file (default file name is tf.mat)"
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input String fileName = "dtf.mat" "Name of the transfer function data file"
        annotation (
          Dialog(
            loadSelector(
              filter="MAT files (*.mat);; All files (*.*)",
              caption="Transfer function data file")));
      input String numName = "n" "Name of the numenator of the transfer function";
      input String denName = "d"
        "Name of the denominator of the transfer function";

    protected
      Integer numSize[2] = Streams.readMatrixSize(fileName, numName) annotation(__Dymola_allowForSize=true);
      Integer denSize[2] = Streams.readMatrixSize(fileName, denName) annotation(__Dymola_allowForSize=true);

      Real num[numSize[1],numSize[2]]=
        Streams.readRealMatrix(fileName, numName, numSize[1], numSize[2])
        "Numenator coefficients";
      Real den[denSize[1],denSize[2]]=
        Streams.readRealMatrix(fileName, denName, denSize[1], denSize[2])
        "Denominator coefficients";
      Integer ns2 = numSize[2] annotation(__Dymola_allowForSize=true);
      Integer ds2 = denSize[2] annotation(__Dymola_allowForSize=true);
      Real Ts[1,1] = Streams.readRealMatrix(fileName, "Ts", 1, 1);

    public
      output DiscreteTransferFunction dtf(n=fill(0,ns2), d=fill(0,ds2))
        "Discrete transfer function";

    algorithm
      dtf.n := vector(num);
      dtf.d := vector(den);
      dtf.uName := numName;
      dtf.yName := denName;
      dtf.Ts := scalar(Ts);

        annotation (Documentation(info="<html>
</html>"));
    end fromFile;

  function fromModel
    "Generate a DiscreteTransferFunction record array from a state space representation resulted from linearization of a model"

    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Simulator = DymolaCommands.SimulatorAPI;

    input String modelName "Name of the Modelica model"  annotation(Dialog(__Dymola_translatedModel(translate=true)));
    input Real T_linearize = 0
        "Point in time of simulation to linearize the model";
    input String fileName = "dslin" "Name of the result file";
    input Modelica.SIunits.Time Ts = 1 "Sample time";
    input Modelica_LinearSystems2.Utilities.Types.Method method=
      Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

    protected
    String fileName2 = fileName + ".mat";
    Boolean OK1 = Simulator.simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
    Boolean OK2 = Simulator.importInitial("dsfinal.txt");
    Boolean OK3 = Simulator.linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize + 1);
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
      redeclare Real D[ny,nu]) "= model linearized at initial point";
    DiscreteStateSpace dss(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu],
      redeclare Real B2[nx,nu]) "= model linearized at initial point";
    DiscreteStateSpace dss_siso(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,1],
      redeclare Real C[1,nx],
      redeclare Real D[1,1],
      redeclare Real B2[nx,1]) "= model linearized at initial point";

    public
    output DiscreteTransferFunction dtf[:,:];

  algorithm
    ss.A := ABCD[1:nx, 1:nx];
    ss.B := ABCD[1:nx, nx + 1:nx + nu];
    ss.C := ABCD[nx + 1:nx + ny, 1:nx];
    ss.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    ss.uNames := xuyName[nx + 1:nx + nu];
    ss.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    ss.xNames := xuyName[1:nx];

    dss := DiscreteStateSpace(ss, Ts=Ts, method=method);
    dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunctionMIMO(dss);

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
    //         dtf[ic, ib] := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss_siso);
    //     end for;
    //   end for;

    annotation (
      Documentation(info="<html>
</html>"));
  end fromModel;

  end Import;

  annotation (Icon(graphics={
        Rectangle(
          lineColor={160,160,164},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid,
          extent={{-100,-100},{100,100}},
          radius=25.0),
        Text(
          lineColor={255,255,170},
          extent={{-90,-50},{90,50}},
          textString="tf")}));
end DiscreteTransferFunction;
