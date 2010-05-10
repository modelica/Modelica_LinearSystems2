within Modelica_LinearSystems2;
record DiscreteTransferFunction
  "Discrete transfer function description of a single input, single output system (data + operations)"

  extends Modelica.Icons.Record;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));

  Modelica.SIunits.Time Ts "Sample time" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method" 
        annotation(Dialog(group="Data used to construct discrete from continuous system"));

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
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n={r}, d={1});

    algorithm
      dtf.Ts := Ts;
      dtf.method := method;
      dtf.uName := uName;
      dtf.yName := yName;
      annotation (overloadsConstructor=true);
    end fromReal;

    encapsulated function fromZerosAndPoles
      "Generate a discrete transfer function from a set of zeros and poles"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.Math.Complex;

      input Complex z[:]=fill(Modelica_LinearSystems2.Math.Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
      input Complex p[:]=fill(Modelica_LinearSystems2.Math.Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
      input Real k=1.0 "Constant multiplied with transfer function";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(redeclare Real n[size(z, 1)+1], redeclare
          Real d[                                                                          size(p, 1)+1])
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
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2;

           input Real n[:] "Coefficients of numerator polynomial";
           input Real d[:] "Coefficients of denominator polynomial";
           input Modelica.SIunits.Time Ts "Sample time";
           input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";

           input String uName = "" "input name";
           input String yName = "" "output name";

           output DiscreteTransferFunction dtf(redeclare Real n[size(n, 1)], redeclare
          Real d[                                                                             size(d, 1)])
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
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial n "Numerator polynomial";
      input Polynomial d "Denominator polynomial";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n=n.c, d=d.c, Ts=Ts, method=method, uName=uName, yName=yName);

    algorithm
      annotation (overloadsConstructor=true);
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

      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method" 
            annotation(Dialog(group="Data used to construct discrete from continuous system"));

      output DiscreteTransferFunction dtf;
    protected
      StateSpace ss = StateSpace(tf);
      Modelica_LinearSystems2.DiscreteStateSpace dss=
                               Modelica_LinearSystems2.DiscreteStateSpace(ss,Ts,method);

    algorithm
      dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss);
      annotation (overloadsConstructor=true);
    end fromTransferFunction;
  end 'constructor';

 encapsulated operator '-'
    function subtract "Subtract two discrete transfer functions (dtf1 - dtf2)"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf1;
      input DiscreteTransferFunction dtf2;

      output DiscreteTransferFunction result;

    protected
      TransferFunction tf=(Polynomial(dtf1.n)*Polynomial(dtf2.d) - Polynomial(dtf2.n)
          *Polynomial(dtf1.d))/(Polynomial(dtf1.d)*Polynomial(dtf2.d));

    algorithm
      assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"-.subtract\".");
      result := DiscreteTransferFunction(tf, dtf1.Ts, dtf1.method);
    end subtract;

    function negate "Unary minus (multiply transfer function by -1)"
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input DiscreteTransferFunction dtf;
      output DiscreteTransferFunction result(n=-dtf.n, d=dtf.d, Ts=dtf.Ts, method=dtf.method) "= -dtf";
    algorithm
    end negate;
 end '-';

  encapsulated operator function '+'
    "Parallel connection of two discrete transfer functions (= inputs are the same, outputs of the two systems are added)"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

    input DiscreteTransferFunction dtf1 "Transfer function system 1";
    input DiscreteTransferFunction dtf2 "Transfer function system 2";
    output DiscreteTransferFunction result;

  protected
    TransferFunction tf=(Polynomial(dtf1.n)*Polynomial(dtf2.d) + Polynomial(
        dtf2.n)*Polynomial(dtf1.d))/(Polynomial(dtf1.d)*Polynomial(dtf2.d));
  algorithm
    assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps,
      "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"+\".");
    result := DiscreteTransferFunction(tf, Ts=dtf1.Ts, method=dtf1.method);
  end '+';

  encapsulated operator function '*'
    "Multiply two DiscreteTransferFunctions (dtf1 * dtf2)"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

    input DiscreteTransferFunction dtf1 "Transfer function system 1";
    input DiscreteTransferFunction dtf2 "Transfer function system 1";
    output DiscreteTransferFunction result;

  protected
    TransferFunction tf=(Polynomial(dtf1.n)*Polynomial(dtf2.n))/(Polynomial(dtf1.d)
        *Polynomial(dtf2.d));
  algorithm
    assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"+\".");
    result := DiscreteTransferFunction(tf, Ts=dtf1.Ts, method=dtf1.method);
  end '*';

  encapsulated operator function '/'
    "Divide two transfer functions (dtf1 / dtf2)"
     import Modelica;
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.TransferFunction;
     import Modelica_LinearSystems2.DiscreteTransferFunction;

     input DiscreteTransferFunction dtf1 "Transfer function system 1";
     input DiscreteTransferFunction dtf2 "Transfer function system 1";
     output DiscreteTransferFunction result;
  protected
    TransferFunction tf=(Polynomial(dtf1.n)*Polynomial(dtf2.d))/(Polynomial(dtf1.d)*
      Polynomial(dtf2.n));

  algorithm
    assert(abs(dtf1.Ts - dtf2.Ts) <= Modelica.Constants.eps, "Two discrete transfer function systems must have the same sample time Ts for subtraction with \"/\".");
    result := DiscreteTransferFunction(tf,Ts=dtf1.Ts, method=dtf1.method);

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
    TransferFunction tf = (Polynomial(dtf.n)^k)/(Polynomial(dtf.d)^k);

  algorithm
    result := DiscreteTransferFunction(tf, Ts=dtf.Ts, method=dtf.method);
  end '^';

  encapsulated operator function '=='
    "Check whether two transfer functions are identical"

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
    "Transform TransferFunction into a String representation"
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
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

  output DiscreteTransferFunction dtf(n={1,0}, d={1}) "z";
algorithm

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right> z </td><td align=center> =  </td>  <td> DiscreteTransferFunction.<b>s</b>()  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Generate the complex variable z=exp(T*s) as a DiscreteTransferFunction. It can be used for generating like 
<blockquote><pre>
        DiscreteTransferFunction dtf = z/(3*z^2 + 2*z +2)
</pre></blockquote>

</p>
 

 
 
</html> "));
end z;

encapsulated package Analysis

encapsulated function timeResponse
      "Calculate the time response of a zeros-and-poles transfer function"

   import Modelica;
   import Modelica_LinearSystems2;
   import Modelica_LinearSystems2.DiscreteStateSpace;
   import Modelica_LinearSystems2.DiscreteTransferFunction;
   import Modelica_LinearSystems2.Types.TimeResponse;

 extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;     // Input/Output declarations of discrete time response functions
 input Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

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

  (y,t,x_discrete) := DiscreteStateSpace.Analysis.timeResponse(dss=dss, tSpan=tSpanVar, response=response, x0=x0);

   annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> TransferFunction.Analysis.<b>timeResponse</b>(tf, dt, tSpan, responseType, x0)  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
First, the transfer function representation is transformed into state space representation which is given to StateSpace.Analysis.timeResponse and the time response of the state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. 
<blockquote><pre>
    Impulse \"Impulse response\",
    Step \"Step response\",
    Ramp \"Ramp response\",
    Initial \"Initial condition response\"
</pre></blockquote>
The state space system is transformed to a appropriate discrete state space system and, starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and x are calculated for each time step t=k*dt.
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=1/(s^2+s+1);

  Real Ts=0.1;
  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[1]={0,0};
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 
 
<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(tf,Ts,tSpan,response,x0);
//  y[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
</pre></blockquote>
 
</html>"));
end timeResponse;

encapsulated function impulseResponse "Calculate the impulse time response"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

algorithm
  (y,t,x_discrete) := Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(
      dtf=dtf,
      tSpan=0,
      response=Modelica_LinearSystems2.Types.TimeResponse.Impulse,
      x0=zeros(Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> TransferFunction.Analysis.<b>timeResponse</b>(tf, dt, tSpan, responseType, x0)  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
First, the transfer function representation is transformed into state space representation which is given to StateSpace.Analysis.timeResponse and the time response of the state space system is calculated. The type of the time response is defined by the input <b>responseType</b>, i.e. 
<blockquote><pre>
    Impulse \"Impulse response\",
    Step \"Step response\",
    Ramp \"Ramp response\",
    Initial \"Initial condition response\"
</pre></blockquote>
The state space system is transformed to a appropriate discrete state space system and, starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and x are calculated for each time step t=k*dt.
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=1/(s^2+s+1);

  Real Ts=0.1;
  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[1]={0,0};
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 
 
<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(tf,Ts,tSpan,response,x0);
//  y[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0048, 0.0187, 0.04, 0.0694}
</pre></blockquote>
 
 
</html> "));
end impulseResponse;

encapsulated function stepResponse "Calculate the step time response"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

algorithm
  (y,t,x_discrete) := Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(
      dtf=dtf,
      tSpan=0,
      response=Modelica_LinearSystems2.Types.TimeResponse.Step,
      x0=zeros(Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> TransferFunction.Analysis.<b>stepResponse</b>(tf, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>stepResponse</b> calculates the step response of a transfer function. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
TransferFunction.Analysis.stepResponse(tf, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Step, x0=fill(0,TransferFunction.Analysis.denominatorDegree(tf))).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=1/(s^2+s+1);

  Real Ts=0.1;
  Real tSpan= 0.4;
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 
 
<b>algorithm</b>
  (y,t,x):=TransferFunction.Analysis.stepResponse(tf,Ts,tSpan);
//  y[:,1,1]={0, 0.0048, 0.01867, 0.04, 0.0694}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0048, 0.01867, 0.04, 0.0694}
</pre></blockquote>
  
</html> "));
end stepResponse;

encapsulated function rampResponse "Calculate the ramp time response"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

algorithm
  (y,t,x_discrete) := Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(
      dtf=dtf,
      tSpan=0,
      response=Modelica_LinearSystems2.Types.TimeResponse.Ramp,
      x0=zeros(Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.denominatorDegree(dtf)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> TransferFunction.Analysis.<b>rampResponse</b>(ss, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>rampResponse</b> calculates the time response of a transfer function for ramp imput u = t. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
TransferFunction.Analysis.rampResponse(ss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,TransferFunction.Analysis.denominatorDegree(tf))).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=1/(s^2+s+1);

  Real Ts=0.1;
  Real tSpan= 0.4;
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 
 
<b>algorithm</b>
  (y,t,x):=TransferFunction.Analysis.rampResponse(tf,Ts,tSpan);
//  y[:,1,1]={0, 0.0002, 0.0012, 0.0042, 0.0096}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0002, 0.0012, 0.0042, 0.0096}
</pre></blockquote>
 
 
 
</html> "));
end rampResponse;

encapsulated function initialResponse "Calculate the initial time response"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

  input Real x0[:]=fill(0,0) "Initial state vector";

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_tf_discrete;

algorithm
  (y,t,x_discrete) := Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.timeResponse(
      dtf=dtf,
      tSpan=0,
      response=Modelica_LinearSystems2.Types.TimeResponse.Initial,
      x0=x0);

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> TransferFunction.Analysis.<b>initialResponse</b>(tf, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initialResponse</b> calculates the time response of a state space system for given initial condition and zero inputs. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
TransferFunction.Analysis.initialResponse(x0,tf, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=1/(s^2+s+1);

  Real Ts=0.1;
  Real tSpan= 0.4;
  Real x0[2] = {1,1};
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 
 
<b>algorithm</b>
  (y,t,x):=TransferFunction.Analysis.initialResponse(x0,tf,Ts,tSpan);
//  y[:,1,1]={1, 1.0903, 1.1616, 1.2151, 1.252}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 1.0903, 1.1616, 1.2151, 1.252}
</pre></blockquote>
 
 
 
</html> "));
end initialResponse;

encapsulated function denominatorDegree "Return denominator degree"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

  input DiscreteTransferFunction dtf "discrete transfer function of a system";
  output Integer result;

algorithm
  result := size(dtf.d,1)-1;
  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Analysis.<b>denominatorDegree</b>(dtf)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of a discrete transfer function. 

</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();
   Modelica_LinearSystems2.DiscreteTransferFunction dtf=(z+1)/(z^2+z+1);
 
   Real dDegree;

<b>algorithm</b>
  dDegree := TransferFunction.Analysis.denominatorDegree(dtf);
//  dDegree = 2
</pre></blockquote>


</html> "));
end denominatorDegree;

end Analysis;

encapsulated package Conversion

encapsulated function toDiscreteZerosAndPoles
      "Generate a DiscreteZerosAndPoles object from a DiscreteTransferFunction object"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles.Internal;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Math.Complex;

  input DiscreteTransferFunction dtf "transfer function of a system";
  output Modelica_LinearSystems2.DiscreteZerosAndPoles dzp(
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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dzp </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Conversion.<b>toDiscreteZerosAndPoles</b>(tf)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a DiscreteZerosAndPoles record
 <blockquote><pre>
                 product(z + n1[i]) * product(z^2 + n2[i,1]*z + n2[i,2])
        zp = k*---------------------------------------------------------
                product(z + d1[i]) * product(z^2 + d2[i,1]*z + d2[i,2])
</pre></blockquote>of a discrete transfer function representated by numerator and denominator polynomial. The poles and zeros and the gain <tt>k</tt> are computed by
(<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs the DiscreteZerosAndPoles constructor.


<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   DiscreteTransferFunction z = Modelica_LinearSystems2.DiscreteTransferFunction.z();  
   Modelica_LinearSystems2.DiscreteTransferFunction dtf = 1/(z^2 + 3*z +2)


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
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteStateSpace;
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
   dss.C := {b[1:nx] - d*a[1:nx]};

end if;
  dss.D := [d];
  dss.Ts := dtf.Ts;
  dss.method := dtf.method;

 annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> DiscreteTransferFunction.Conversion.toStateSpace<b>toDiscreteStateSpace</b>(dtf)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Transforms a discrete transfer function into discrete state space representation.
There are an infinite number of possible realizations.
Here, the transfer function is transformed into
controller canonical form, i.e. the transfer function
<blockquote><pre>
       b4*z^4 + b3*z^3 + b2*z^2 + b1*z + b0
  y = -------------------------------------- *u
       a4*z^4 + a3*z^3 + a2*z^2 + a1*z + a0
</pre></blockquote>
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
If the numerator polynomial is 1, then the state vector
<b>x</b> is built up of the y(k) (the privious y) and of all the nx-1 predecessor
(nx is the dimension of the state vector):
<blockquote><pre>
   <b>x</b>(k+1) = {y(k-n+1), y(k-n+2), ..., y(k)};
</pre></blockquote>
Note, the state vector <b>x</b> of Modelica.Blocks.Continuous.TransferFunction
is defined slightly differently.

</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
end toDiscreteStateSpace;

end Conversion;

encapsulated package Plot "Functions to plot state space system responses"

encapsulated function bode "Plot transfer function as bode plot"
      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Complex;

      import Modelica_LinearSystems2.Utilities.Plot;

      import SI = Modelica.SIunits;

  input DiscreteTransferFunction dtf "DiscreteTransfer function to be plotted";
  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false" 
                                                    annotation(Dialog(enable=not autoRange));
  input Modelica.SIunits.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                              annotation(Dialog(enable=not autoRange));

  input Boolean magnitude=true "= true, to plot the magnitude of tf" 
                                                                    annotation(choices(__Dymola_checkBox=true));
  input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(__Dymola_checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot of  dtf = "
         + String(dtf)));

    protected
  SI.AngularVelocity w[nPoints];
  Complex z[nPoints];
  SI.Frequency f[nPoints];
  Modelica.SIunits.Conversions.NonSIunits.Angle_deg phi[nPoints];
  Real A[nPoints];
  Boolean OK;
  Complex c;
  Modelica.SIunits.Angle phi_old;
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
    denZeros[i] := Complex.log(denZerosZ[i])/dtf.Ts;
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
    w[i] := Modelica.SIunits.Conversions.from_Hz(f[i]);
    z[i] := Complex.exp(Complex(0,w[i]*dtf.Ts));
    c := TransferFunction.Analysis.evaluate(
          tf,
          z[i],
          1e-10);
    A[i] := Complex.'abs'(c);
    phi_old := Complex.arg(c, phi_old);
    phi[i] := Modelica.SIunits.Conversions.to_deg(phi_old);

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

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>plotBode</b>(tf)
   or
TransferFunction.Plot.<b>plotBode</b>(tf, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Plots the bode-diagram of a transfer function.


</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf =(s^2 + 5*s + 7)/(s^2 + 5*s + 6);
   
<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.plotBode(tf)
//  gives:
</pre></blockquote>

</p>
 
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhase.png\">
 
</blockquote>
<p>


</html> "));
end bode;

encapsulated function timeResponse
      "Plot the time response of a system represented by a transfer function. The response type is selectable"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

  input Modelica_LinearSystems2.DiscreteTransferFunction dtf;
//  input Real dt=0 "Sample time [s]";
  input Real tSpan=0 "Simulation time span [s]";

  input Modelica_LinearSystems2.Types.TimeResponse response=
      Modelica_LinearSystems2.Types.TimeResponse.Step "Type of time response";
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
  annotation (interactive=true, Documentation(info="<html>
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
TransferFunction.Plot.<b>timeResponse</b>(tf);
   or
TransferFunction.Plot.<b>timeResponse</b>(tf, dt, tSpan,response, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
                   device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>


<p><b><font style=\"color: #008000; \">Description</font></b></p>
<p>Function <b>timeResponse</b> plots the time response of a transfer function. The character of the time response if defined by the input 
<a href=\"Modelica://Modelica_LinearSystems2.Types.TimeResponse\">response</a>, i.e. Impulse, Step, Ramp, or Initial. See also 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>, and 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>. </p>

<p><b><font style=\"color: #008000; \">Example</font></b></p>
<pre>   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  </pre>
<pre>   Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);</pre>
<pre><br/>   Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;</pre>
<pre><br/><b>algorithm</b></pre>
<pre>   Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(tf, dt=0.02, tSpan=3, response=response)</pre>
<pre>//  gives: </pre>
<p><img src=\"modelica://Modelica_LinearSystems2/Extras/Images/timeResponse.png\"/> </p>
</html>"));
end timeResponse;

encapsulated function impulse "Impulse response plot"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

    import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteTransferFunction dtf "zeros-and-poles transfer function";
    input Real tSpan=0 "Simulation time span [s]";

    input Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
        DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
         Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response of  zp = "
           + String(dtf)));

    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Impulse
        "type of time response";

algorithm
// set sample time
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.timeResponse(
      dtf=dtf,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>impulse</b>(tf)  
   or
TransferFunction.Plot.<b>impulse</b>(tf, dt, tSpan, x0, columnLabels, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>impulse</b> plots the impulse response of a transfer function. It is based on <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a> . See also
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>

   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);

<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.impulse(tf, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/impulseResponse.png\">
</p>
<p>
</p>

</html> "));
end impulse;

encapsulated function step "Step response plot"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteTransferFunction dtf;
  input Real tSpan=0 "Simulation time span [s]";

  input Modelica_LinearSystems2.Types.TimeResponse response=
      Modelica_LinearSystems2.Types.TimeResponse.Step "type of time response";
  input Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
      DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  dtf = "
         + String(dtf)));

algorithm
 DiscreteTransferFunction.Plot.timeResponse(
    dtf=dtf,
    tSpan=tSpan,
    response=response,
    x0=x0,
    defaultDiagram=defaultDiagram,
    device=device);

equation

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>step</b>(tf)  
   or
TransferFunction.Plot.<b>step</b>(tf, dt, tSpan, x0, columnLabels, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>step</b> plots the step response of a transfer function. It is based on <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a> . See also
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>

   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);

<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.step(tf, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/stepResponse.png\">
</p>
<p>
</p>

</html>"));
end step;

encapsulated function ramp "Ramp response plot"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteTransferFunction dtf;
  input Real tSpan=0 "Simulation time span [s]";

  input Modelica_LinearSystems2.Types.TimeResponse response=
      Modelica_LinearSystems2.Types.TimeResponse.Ramp "type of time response";
  input Real x0[DiscreteTransferFunction.Analysis.denominatorDegree(dtf)]=zeros(
      DiscreteTransferFunction.Analysis.denominatorDegree(dtf))
        "Initial state vector";

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  dtf = "
         + String(dtf)));

algorithm
 Modelica_LinearSystems2.DiscreteTransferFunction.Plot.timeResponse(
    dtf=dtf,
    tSpan=tSpan,
    response=response,
    x0=x0,
    defaultDiagram=defaultDiagram,
    device=device);

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>ramp</b>(tf)  
   or
TransferFunction.Plot.<b>ramp</b>(tf, dt, tSpan, x0, columnLabels, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>ramp</b> plots the ramp response of a transfer function. It is based on <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a> . See also
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>.


</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf =(2*s^2 + 7*s + 13)/(s^3 + 6*s^2 + 17*s + 12);

<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.ramp(tf)
//  gives:
</pre></blockquote>

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/rampResponse.png\">
</p>
<p>
</p>

</html> "));
end ramp;

encapsulated function initialResponse "Initial condition response plot"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

  input Modelica_LinearSystems2.DiscreteTransferFunction dtf;
  input Real tSpan=0 "Simulation time span [s]";

  input Modelica_LinearSystems2.Types.TimeResponse response=
      Modelica_LinearSystems2.Types.TimeResponse.Initial
        "type of time response";
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

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>initialResponse</b>(tf)  
   or
TransferFunction.Plot.<b>initialResponse</b>(tf, dt, tSpan, y0, columnLabels, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initialResponse</b> plots the initial response, i.e. the zeros input response of a transfer function. It is based on <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a> . See also
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">initialResponse</a>.




</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf = (s + 1)/(s^2 + 5*s + 12);
   Real y0=1; 



<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.initialResponse(tf,y0=y0, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/initialResponseTF.png\">
</p>
<p>
</p>

</html> "));
end initialResponse;

end Plot;

encapsulated package Import

function fromModel
      "Generate a TransferFunction record array from a state space representation resulted from linearization of a model"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.DiscreteTransferFunction;

  input String modelName "Name of the Modelica model" annotation(Dialog(translatedModel));
  input Real T_linearize=0 "point in time of simulation to linearize the model";
  input String fileName="dslin" "Name of the result file";
  input Modelica.SIunits.Time Ts=1 "Sample time";
  input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";

    protected
  String fileName2=fileName + ".mat";
  Boolean OK1=simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
  Boolean OK2=importInitial("dsfinal.txt");
  Boolean OK3=linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize + 1);
  Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
  Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
  Integer nx=integer(nxMat[1, 1]);
  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
  String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);

  StateSpace result(
    redeclare Real A[nx,nx],
    redeclare Real B[nx,nu],
    redeclare Real C[ny,nx],
    redeclare Real D[ny,nu]) "= model linearized at initial point";
  TransferFunction tf[:,:];
    public
  output DiscreteTransferFunction dtf[:,:];

algorithm
  result.A := ABCD[1:nx, 1:nx];
  result.B := ABCD[1:nx, nx + 1:nx + nu];
  result.C := ABCD[nx + 1:nx + ny, 1:nx];
  result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
  result.uNames := xuyName[nx + 1:nx + nu];
  result.yNames := xuyName[nx + nu + 1:nx + nu + ny];
  result.xNames := xuyName[1:nx];

  tf := Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunctionMIMO(
    result);

    for i in 1:ny loop
      for j in 1:nu loop
        dtf[i,j] := DiscreteTransferFunction(tf=tf[i,j], Ts=Ts, method=method);
      end for;
    end for;

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> TransferFunction.Import.<b>fromModel</b>(modelName, T_linearize, fileName)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Generate a matrix of TransferFunction data records by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation. The system is genrated by using <a href=\"Modelica://Modelica_LinearSystems2.
StateSpace.Import.fromFile\">StateSpace.Import.fromFile</a> followed by a conversion from sate space to transfer function representation.
 
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   String modelName = \"Modelica_LinearSystems2.Examples.DoublePendulum\"; 
   Real T_linearize = 5;; 
   
 
<b>algorithm</b>
  tf = Modelica_LinearSystems2.TransferFunction.Import.fromModel(modelName, T_linearize);
 
//  tf = [(0.13*s^4 + 0.05558*s^3 + 1.12241*s^2 - 5.16971*s + 9.04744)/(s^6 + 0.09*s^5 + 9.13717*s^4 - 32.0637*s^3 + 58.78*s^2 + 6.3659e-014*s - 1.1703e-014);
          (0.13*s^4 + 0.05558*s^3 + 1.12241*s^2 - 5.16971*s + 9.04744)/(s^5 + 0.09*s^4 + 9.13717*s^3 - 32.0637*s^2 + 58.78*s - 2.7929e-015);
          (-0.014*s^2 + 0.31906*s - 0.8106)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.014*s^3 + 0.31906*s^2 - 0.8106*s)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.1*s^2 - 0.160918*s - 0.21842)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.1*s^3 - 0.160918*s^2 - 0.21842*s)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78)]
                      
</pre></blockquote>
 
 
 
 
</html> "));
end fromModel;

encapsulated function fromFile
      "Generate a transfer function data record by reading numenator coefficients and denominator coefficients from a file (default file name is tf.mat)"

    import Modelica_LinearSystems2.DiscreteTransferFunction;
    import Modelica_LinearSystems2.Math.Polynomial;
  input String fileName="tf.mat" "Name of the transfer function data file"   annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="transfer function data file")));
  input String numName="n" "Name of the numenator of the transfer function";
  input String denName="d" "Name of the denominator of the transfer function";

    protected
  input Integer numSize[2]=readMatrixSize(fileName, numName);
  input Integer denSize[2]=readMatrixSize(fileName, denName);

  Real num[numSize[1],numSize[2]]=readMatrix(
        fileName,
        numName,
        numSize[1],
        numSize[2]) "numenator coefficients";
  Real den[denSize[1],denSize[2]]=readMatrix(
        fileName,
        denName,
        denSize[1],
        denSize[2]) "denominator coefficients";
  input Integer ns2=numSize[2];
  input Integer ds2=denSize[2];
  Real Ts[1,1]=readMatrix(fileName, "Ts", 1, 1);

    public
 output DiscreteTransferFunction dtf(n=fill(0,ns2),d=fill(0,ds2))
        "Discrete transfer function";

algorithm
  dtf.n := vector(num);
  dtf.d := vector(den);
  dtf.uName := numName;
  dtf.yName := denName;
  dtf.Ts := scalar(Ts);

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> TransferFunction.Import.<b>fromFile</b>(fileName, numName, denName)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Reads and loads a transfer function from a mat-file <tt>fileName</tt>. The file must contain the names of the vector with the polynomial coefficients of numerator and denominator

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
     

<b>algorithm</b>
  tf:=Modelica_LinearSystems2.TransferFunction.Import.fromFile(\"tf.mat\", \"n\", \"d\");
//  tf = (s^2 + 2*s + 3)/(4*s^2 + 5*s + 6)
</pre></blockquote>


</html> "));
end fromFile;

end Import;

end DiscreteTransferFunction;
