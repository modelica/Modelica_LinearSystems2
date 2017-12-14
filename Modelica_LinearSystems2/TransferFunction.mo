within Modelica_LinearSystems2;
operator record TransferFunction
  "Continuous transfer function description of a single input, single output system (data + operations)"

  import Modelica_LinearSystems2.Math.Polynomial;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));

  String uName="u" "Name of input signal" annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal" annotation(Dialog(group="Signal names"));

/* If the numerator polynomial has no coefficients, the transfer function
   is zero. The denominator polynomial must always have at
   least one coefficient, such as {1}
*/

  encapsulated operator 'constructor'
    "Collection of operators to construct a TransferFunction data record"
    import Modelica;
    import Modelica_LinearSystems2.TransferFunction;

    function fromReal
      "Generate a TransferFunction data record from a real value"
      import Modelica;
      import Modelica_LinearSystems2.TransferFunction;

      input Real r "Value of Real variable";
      input String uName="" "input name";
      input String yName="" "output name";
      output TransferFunction tf(n={r}, d={1});

    algorithm
      tf.uName := uName;
      tf.yName := yName;
    end fromReal;

  encapsulated function fromZerosAndPoles
    "Generate a TransferFunction data record from a set of zeros and poles"

    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Complex;

    input Complex z[:]=fill(Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
    input Complex p[:]=fill(Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
    input Real k=1.0 "Constant multiplied with transfer function";
    input String uName="" "input name";
    input String yName="" "output name";
    output TransferFunction tf(redeclare Real n[size(z, 1)+1], redeclare Real d[size(p, 1)+1])
        "TransferFunction built by ZerosAndPoles object";

    protected
    Polynomial pn = k*Polynomial(z);
    Polynomial pd = Polynomial(p);

  algorithm
    tf.n := pn.c;
    tf.d := pd.c;
    tf.uName := uName;
    tf.yName := yName;

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
      "Generate a TransferFunction data record from numerator and denominator array"
      import Modelica;
      import Modelica_LinearSystems2.TransferFunction;

           input Real n[:] "Coefficients of numerator polynomial";
           input Real d[:] "Coefficients of denominator polynomial";
           input String uName = "" "input name";
           input String yName = "" "output name";

           output TransferFunction tf(redeclare Real n[size(n, 1)], redeclare Real
               d[                                                                    size(d, 1)])
        "Transfer function";

    algorithm
                 //this is the constructor algorithm
           assert(size(d, 1) > 0, "Input denominator d must have at least one element, however\n"
              + "d is an empty vector");
           tf.n := n;
           tf.d := d;
           tf.uName := uName;
           tf.yName := yName;
    end fromArrays;

    function fromPolynomials
      "Generate a TransferFunction data record from a numerator and denominator polynomial"
      import Modelica;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial n "Numerator polynomial";
      input Polynomial d "Denominator polynomial";
      input String uName="" "input name";
      input String yName="" "output name";
      output TransferFunction tf(n=n.c, d=d.c,uName=uName, yName=yName);

    algorithm
    end fromPolynomials;

  end 'constructor';

  encapsulated operator '-'
    "Collection of operators for subtraction of transfer functions"

    import Modelica;
    extends Modelica.Icons.Package;
    function subtract "Subtract two transfer functions (tf1 - tf2)"
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf1;
      input TransferFunction tf2;

      output TransferFunction result = (Polynomial(tf1.n)*Polynomial(tf2.d) - Polynomial(tf2.n)*Polynomial(tf1.d))/(Polynomial(tf1.d)*Polynomial(tf2.d));

    algorithm
    end subtract;

    function negate "Unary minus (multiply transfer function by -1)"
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf;
      output TransferFunction result(n=-tf.n, d=tf.d) "= -tf";
    algorithm
    end negate;
    annotation (Documentation(info="<html>
<p>This package contains operators for subtraction of transfer function records. </p>
</html>"));
  end '-';

  encapsulated operator function '+'
    "Parallel connection of two transfer functions (= inputs are the same, outputs of the two systems are added)"
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.TransferFunction;

     input TransferFunction tf1 "Transfer function system 1";
     input TransferFunction tf2 "Transfer function system 2";
     output TransferFunction result;
  algorithm
      result := (Polynomial(tf1.n)*Polynomial(tf2.d) + Polynomial(tf2.n)*Polynomial(tf1.d))/(Polynomial(tf1.d)*Polynomial(tf2.d));
  end '+';

  encapsulated operator function '*'
    "Multiply two TransferFunctions (tf1 * tf2)"
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.TransferFunction;

     input TransferFunction tf1 "Transfer function 1";
     input TransferFunction tf2 "Transfer function 2";
     output TransferFunction result;
  algorithm
     result := (Polynomial(tf1.n)*Polynomial(tf2.n))/(Polynomial(tf1.d)*Polynomial(tf2.d));
  end '*';

  encapsulated operator function '/'
    "Divide two transfer functions (tf1 / tf2)"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.TransferFunction;

    input TransferFunction tf1 "Transfer function system 1";
    input TransferFunction tf2 "Transfer function system 2";
    output TransferFunction result "Result = tf1/tf2";

  algorithm
    result := (Polynomial(tf1.n)*Polynomial(tf2.d))/(Polynomial(tf1.d)*
      Polynomial(tf2.n));
  end '/';

  encapsulated operator function '^'
    "Integer power of TransferFunction (tf1^k)"
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.TransferFunction;

     input TransferFunction tf "Transfer function";
     input Integer k(min=0) = 1 "Integer exponent";
     output TransferFunction result;
  algorithm
    result := (Polynomial(tf.n)^k)/(Polynomial(tf.d)^k);
  end '^';

  encapsulated operator function '=='
    "Check whether two transfer functions are identical"
     import Modelica_LinearSystems2.Math.Polynomial;
     import Modelica_LinearSystems2.TransferFunction;

     input TransferFunction tf1 "Transfer function system 1";
     input TransferFunction tf2 "Transfer function system 1";
     input Real eps(min=0) = 0
      "Two coefficients c1 and c2 of the two transfer functions are identical if abs(c1-c2) <= eps";

     output Boolean result "= tf1 == tf2";
  algorithm
    result := (Polynomial(tf1.n) == Polynomial(tf2.n)) and (Polynomial(tf1.d) == Polynomial(tf2.d));
  end '==';

  encapsulated operator function 'String'
    "Transform TransferFunction into a String representation"
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf
      "Transfer function to be transformed in a String representation";
      input Integer significantDigits=6
      "Number of significant digits that are shown";
      input String name="s" "Independent variable name used for printing";
      output String s="";
  protected
      Integer n_num=size(tf.n, 1) - 1;
      Integer n_den=size(tf.d, 1) - 1;
      Boolean numParenthesis;
  algorithm
      if n_num == -1 then
        s := "0";
      else
        numParenthesis := n_num > 0 and not (n_den == 0 and tf.d[1] == 1);
        if numParenthesis then
          s := "(";
        end if;
         s := s + String(
              Polynomial(tf.n),
              significantDigits,
              name);

        if numParenthesis then
          s := s + ")";
        end if;
      if n_den > 0 or tf.d[1] <> 1 then
          if n_den > 0 then
            s := s + "/(";
          else
            s := s + "/";
          end if;

          s := s + String(
                Polynomial(tf.d),
                significantDigits,
                name);

          if n_den > 0 then
            s := s + ")";
          end if;
        end if;
      end if;
  end 'String';

  encapsulated function s "Generate the transfer function s"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.TransferFunction;

    output TransferFunction tf(n={1,0}, d={1}) "s";
  algorithm

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
s = TransferFunction.<b>s</b>()
</pre></blockquote>

<h4>Description</h4>
<p>
Generate the complex Laplace variable as a TransferFunction.
It can be used for generating like
</p>
<blockquote><pre>
TransferFunction tf = s/(3*s^2 + 2*s +2)
</pre></blockquote>
</html>"));
  end s;

  encapsulated package Analysis
    "Package of functions to analyse transfer function represented by a TransferFunction record"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
    extends Modelica.Icons.Package;
    function analysis
      "Make a system analysis based on the poles and zeros of the system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Internal.AnalyseOptions;
      import Modelica_LinearSystems2.Internal.AnalyseOptions2;
      import Modelica_LinearSystems2.Internal.Eigenvalue;

      input TransferFunction tf(uName="u", yName="y")
        "transfer function of a system";

      input AnalyseOptions2 analyseOptions2=
          Modelica_LinearSystems2.Internal.AnalyseOptions2(printControllability=
           false, printObservability=false);

      input String fileName="eigenvalues.html"
        "Name of html-file that contains eigenvalue table";

      input String systemName="tf"
        "Name of system (used as heading in html file)";
      input String description="" "Description of system (used in html file)";

    protected
      String dummyFileName = "dummy" + fileName;

      input StateSpace ss=StateSpace(tf);

      Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
          AnalyseOptions(
              plotEigenValues=analyseOptions2.plotEigenValues,
              plotInvariantZeros=analyseOptions2.plotInvariantZeros,
              plotStepResponse=analyseOptions2.plotStepResponse,
              plotFrequencyResponse=analyseOptions2.plotFrequencyResponse,
              printSystem=analyseOptions2.printSystem,
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
      assert(TransferFunction.Analysis.denominatorDegree(tf) >=
        TransferFunction.Analysis.numeratorDegree(tf),
        " Denominator polynominal of transfer function in function\"TransferFunction.Analysis.analysis\"has to be of higher or equal order than numerator polynomial");

     Modelica.Utilities.Files.removeFile(fileName);
     Modelica.Utilities.Files.removeFile(dummyFileName);
     if analyseOptions.printSystem and size(ss.A,1) <= 50 then
    Modelica_LinearSystems2.TransferFunction.Analysis.analysis.printSystem(
            tf,
            fileName,
            systemName,
            description);
    Modelica_LinearSystems2.TransferFunction.Analysis.analysis.printSystem(
            tf,
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
    equation

    public
      encapsulated function printSystem
        "Print the state space system in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.TransferFunction;
        import Modelica_LinearSystems2;

        input TransferFunction tf "transfer function to analyze";
        input String fileName="systemAnalysis.html"
          "File on which the transfer fucntion is written in html format";
        input String systemName="Transfer Function" "name of the system";
        input String description="" "Description of system (used in html file)";
        input String format=".3g" "Format of numbers (e.g. \"20.8e\")";
      protected
        String st=String(tf);

      algorithm
        Modelica.Utilities.Files.removeFile(fileName);
        print("<html>\n<body>\n<p><b>System report</b></p>", fileName);
        print("<p> The system " + systemName + " is defined by</p>", fileName);
        print("G(s) = " + st, fileName);
         if description == "" then
          print("</table> ", fileName);
        else
          print("</table>", fileName);
          print("<p><b>Description</b></p>", fileName);
          print(description, fileName);
        end if;
        print("<br></body></html>",fileName);

      end printSystem;

      annotation (__Dymola_interactive=true);
    end analysis;

   encapsulated function timeResponse
      "Calculate the time response of a transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

    extends Modelica_LinearSystems2.Internal.timeResponseMask2_tf;     // Input/Output declarations of time response functions
      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step;

    input Real x0[TransferFunction.Analysis.denominatorDegree(tf)]=zeros(TransferFunction.Analysis.denominatorDegree(tf))
        "Initial state vector";

    protected
     StateSpace ss=StateSpace(tf);

   algorithm
     (y,t,x_continuous) := StateSpace.Analysis.timeResponse(sc=ss, dt=dt, tSpan=tSpan, response=response, x0=x0);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = TransferFunction.Analysis.<b>timeResponse</b>(tf, dt, tSpan, responseType, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the transfer function representation is transformed into state
space representation which is given to StateSpace.Analysis.timeResponse
and the time response of the state space system is calculated. The type
of the time response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step &quot;Step response&quot;,
Ramp &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;.
</pre></blockquote>
<p>
The state space system is transformed to a appropriate discrete state space
system and, starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and
x are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
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
      import Modelica_LinearSystems2.TransferFunction;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_tf;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(
          tf=tf,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
          x0=zeros(Modelica_LinearSystems2.TransferFunction.Analysis.denominatorDegree(tf)));

    annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = TransferFunction.Analysis.<b>timeResponse</b>(tf, dt, tSpan, responseType, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
First, the transfer function representation is transformed into state space
representation which is given to StateSpace.Analysis.timeResponse and the
time response of the state space system is calculated. The type of the time
response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
    Impulse \"Impulse response\",
    Step \"Step response\",
    Ramp \"Ramp response\",
    Initial \"Initial condition response\"
</pre></blockquote>
<p>
The state space system is transformed to a appropriate discrete state space
system and, starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y
and x are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
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
  end impulseResponse;

  encapsulated function stepResponse "Calculate the step time response"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_tf;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(
          tf=tf,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
          x0=zeros(Modelica_LinearSystems2.TransferFunction.Analysis.denominatorDegree(tf)));

    annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = TransferFunction.Analysis.<b>stepResponse</b>(tf, dt, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>stepResponse</b> calculates the step response of a transfer function.
The state space system is transformed to a appropriate discrete state space
system and, starting at
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
TransferFunction.Analysis.stepResponse(tf, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Step, x0=fill(0,TransferFunction.Analysis.denominatorDegree(tf))).
</pre></blockquote>

<h4>Example</h4>
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

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
</html>"));
  end stepResponse;

  encapsulated function rampResponse "Calculate the ramp time response"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_tf;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(
          tf=tf,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
          x0=zeros(Modelica_LinearSystems2.TransferFunction.Analysis.denominatorDegree(tf)));

    annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = TransferFunction.Analysis.<b>rampResponse</b>(tf, dt, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rampResponse</b> calculates the time response of a transfer
function for ramp imput u = t. The state space system is transformed
to a appropriate discrete state space system and, starting at
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
TransferFunction.Analysis.rampResponse(ss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,TransferFunction.Analysis.denominatorDegree(tf))).
</pre></blockquote>

<h4>Example</h4>
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

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
</html>"));
  end rampResponse;

  encapsulated function initialResponse
      "Calculate the time response for given initial condition and zero inputs"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;

    input Real x0[:]=fill(0,0) "Initial state vector";

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask2_tf;

  algorithm
    (y,t,x_continuous) :=Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse(
          tf=tf,
          dt=dt,
          tSpan=tSpan,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
          x0=x0);

    annotation(__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = TransferFunction.Analysis.<b>initialResponse</b>(tf, dt, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the time response of a state space system for
given initial condition and zero inputs. The state space system is transformed
to a appropriate discrete state space system and, starting at
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
TransferFunction.Analysis.initialResponse(x0,tf, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
TransferFunction.Analysis.timeResponse(tf, dt, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>

<h4>Example</h4>
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

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.timeResponse\">TransferFunction.Analysis.timeResponse</a>
</p>
</html>"));
  end initialResponse;

    encapsulated function numeratorDegree "Return numerator degree"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";
      output Integer result;

    algorithm
      result := size(tf.n,1)-1;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = TransferFunction.Analysis.<b>numeratorDegree</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>numeratorDegree</b> calculates the degree of the numerator polynomial of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

   Real nDegree;

<b>algorithm</b>
  nDegree := TransferFunction.Analysis.numeratorDegree(tf);
//  nDegree = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.denominatorDegree\">TransferFunction.Analysis.denominatorDegree</a>
</p>
</html>"));
    end numeratorDegree;

    encapsulated function denominatorDegree "Return denominator degree"
      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";
      output Integer result;

    algorithm
      result := size(tf.d,1)-1;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = TransferFunction.Analysis.<b>denominatorDegree</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

   Real dDegree;

<b>algorithm</b>
  dDegree := TransferFunction.Analysis.denominatorDegree(tf);
//  dDegree = 2
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.numeratorDegree\">TransferFunction.Analysis.numeratorDegree</a>
</p>
</html>"));
    end denominatorDegree;

    encapsulated function evaluate
      "Evaluate a transfer function for a given (Complex) value of s"

      import Modelica;
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "Transfer function of a system";
      input Complex s "Value of s where tf shall be evaluated";
      input Real den_min=0 "Value of |denominator(s)| is limited by den_min";
      output Complex result "= tf(s)";

    protected
      Complex den=Polynomial.evaluateComplex(Polynomial(tf.d), s);
      Real abs_den=Modelica.ComplexMath.'abs'(den);
    algorithm
      den := if abs_den >= den_min then den else -abs_den+0*j;
      result := Polynomial.evaluateComplex(Polynomial(tf.n), s)/den;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = TransferFunction.Analysis.<b>evaluate</b>(tf, s)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates a transfer function at a given (complex) value of s.
The transfer function G(s)=N(s)/D(s) is evaluated by calculating the numerator polynomial N(s) and the denominator polynomial D(s).
</p>

<h4>Example</h4>
<blockquote><pre>
   Complex j = Modelica_LinearSystems2.Math.Complex.j();
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
   Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

   Complex result;

<b>algorithm</b>
  result := Modelica_LinearSystems2.TransferFunction.Analysis.evaluate(tf, j+1);
//  result = 0.538462 - 0.307692j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>
</html>"));
    end evaluate;

    encapsulated function zerosAndPoles
      "Calculate zeros and poles of a transfer function"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";

      output Complex z[:]=Polynomial.roots(Polynomial(tf.n))
        "Zeros (Complex vector of numerator zeros)";
      output Complex p[:]=Polynomial.roots(Polynomial(tf.d))
        "Poles (Complex vector of denominator zeros)";
      output Real k
        "Constant multiplied with transfer function that is factorized with zeros and poles";

    protected
      TransferFunction tf2=TransferFunction(z, p);
      Real r;
      Complex s;
      Complex y1;
      Complex y2;
    algorithm
      // Determine an s-value that is neither a zero nor a pole
      r := 1.0;
      for i in 1:size(z, 1) loop
        r := max(r, abs(z[i].re));
      end for;
      for i in 1:size(p, 1) loop
        r := max(r, abs(p[i].re));
      end for;
      r := 2*r;
      s := Complex(r, 0);

      // Evaluate both tf and tf2 and determine k from the quotient
      y1 := TransferFunction.Analysis.evaluate(tf, s);
      y2 := TransferFunction.Analysis.evaluate(tf2, s);
      k := y1.re/y2.re;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(z,p,k) = TransferFunction.Analysis.<b>zerosAndPoles</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the zeros, poles and gain of transfer function.
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.evaluate\">Analysis.evaluate</a> is used to calculate gain <tt>k</tt>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

public
  output Complex z;
  output Complex p;
  output Real k;

<b>algorithm</b>
  (z,p,k)=Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles(tf);
//  z = {-1}
//  p = {-0.5 + 0.866025j, -0.5 - 0.866025j}
//  k = 1
</pre></blockquote>
</html>"));
    end zerosAndPoles;

    function eigenValues
      "Calculate the eigenvalues of a linear transfer function system and write them in a complex vector"
    //encapsulated function eigenValues
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;

      input TransferFunction tf "transfer function of a system";
      output Complex eigval[:] "eigen values of the system";

    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      assert(TransferFunction.Analysis.denominatorDegree(tf) >
        TransferFunction.Analysis.numeratorDegree(tf),
        " Denominator polynominal of transfer function in function\"TransferFunction.Analysis.eigenValues\"has to be of higher order than numerator polynomial");

      eigval := StateSpace.Analysis.eigenValues(ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
eigenvalues = TransferFunction.Analysis.<b>eigenValues</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvalues of the corresponding state space representation of a transfer function. The output is a complex vector containing the eigenvalues. Note, that the conversion of the transfer function does not result in a minimal state space system. Therefore also unobservable and uncontrollable eigenvalues will be calculated.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

  Complex eigenvalues[2];

<b>algorithm</b>
  eigenvalues = Modelica_LinearSystems2.TransferFunction.Analysis.eigenValues(tf);
// eigenvalues = {-0.5 + j*sqrt(3)/2, -0.5 - j*sqrt(3)/2}
</pre></blockquote>
</html>"));
    end eigenValues;

    encapsulated function eigenVectors
      "Calculate the right eigenvectors of the state space system corresponding to a transfer function and write them columnwise in a matrix. Optionally, the eigenvalues are computed"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica.Math.Matrices.LAPACK;

      input TransferFunction tf "transfer function of a system";
      input Boolean onlyEigenvectors=true;
      output Real eigvec[:,:] "eigen values of the system";
      output Complex eigval[:] "eigen values of the system";
    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      assert(TransferFunction.Analysis.denominatorDegree(tf) >
        TransferFunction.Analysis.numeratorDegree(tf),
        " Denominator polynominal of transfer function in function\"TransferFunction.Analysis.eigenVectors\"has to be of higher order than numerator polynomial");
      (eigvec,eigval) := StateSpace.Analysis.eigenVectors(ss=ss,
        onlyEigenvectors=onlyEigenvectors);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(eigenvectors, eigenvalues) = TransferFunction.Analysis.<b>eigenVectors</b>(tf, onlyEigenvectors)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvectors and optionally (onlyEigenvectors=false) the eigenvalues of the corresponding state space system of a transfer function. The output <tt>eigenvectors</tt> is a matrix with the same dimension as matrix <b>ss.A</b>. Just like in <a href=\"modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>, if the i-th eigenvalue has an imaginary part, then <tt>eigenvectors</tt>[:,i] is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the eigenvector of the i-th eigenvalue.<br>
The eigenvalues are returned as a complex vector <tt>eigenvalues</tt>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(2*s+2)/(s^2+2*s+2);

  Real eigenvectors[2,2];
  Complex eigenvalues[2];

<b>algorithm</b>
  (eigenvectors, eigenvalues) = Modelica_LinearSystems2.TransferFunction.Analysis.eigenVectors(tf, true);
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
      "Compute invariant zeros of linear transfer function"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;

      input TransferFunction tf "transfer function of a system";

      output Complex Zeros[:] "invariant zeros";

    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      Zeros := StateSpace.Analysis.invariantZeros(ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zeros = TransferFunction.Analysis.<b>invariantZeros</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes the invariant zeros of the corresponding state space representation of a transfer function. The output is a complex vector containing the eigenvalues. Note, that the conversion of the transfer function does not result in a minimal state space system. Therefore, also zeros equal to unobservable or uncontrollable eigenvalues will be computed.

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

  Complex zeros[:];

<b>algorithm</b>
  zeros := Modelica_LinearSystems2.TransferFunction.Analysis.invariantZeros(tf);
// zeros = {-1}
</pre></blockquote>
</html>"));
    end invariantZeros;

    encapsulated function dcGain
      "Return steady state gain k (for a stable system: k = value of y at infinite time for a step input)"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "Transfer function of a system";
      output Real k "Steady state gain";
      output Boolean finite = true
        "True, if k is finite, otherwise k is infinite (k=Modelica.Constants.inf returned)";
    protected
      StateSpace ss=StateSpace(tf);
      Real K[1,1];
    algorithm
      (K, finite) := StateSpace.Analysis.dcGain(ss=ss);
      k :=K[1, 1];

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
          k = <b>dcGain</b>(tf);
(k, finite) = <b>dcGain</b>(tf);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the steady state gain <b>k</b> of a
TransferFunction tf(s), i.e. k = tf(s=0).
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
      "Check controllability of a transfer function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

        input TransferFunction tf "transfer function of a system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

        output Boolean controllable;
    protected
        StateSpace ss=StateSpace(tf);

    algorithm
        controllable := StateSpace.Analysis.isControllable(ss=ss, method=method);

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
controllable = TransferFunction.Analysis.<b>isControllable</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function TransferFunction.Analysis.<b>isControllable</b> checks the controllability of a transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isControllable\">StateSpace.Analysis.isControllable</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2 + 2*s +1);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean controllable;

<b>algorithm</b>
  controllable := Modelica_LinearSystems2.StateSpace.Analysis.isControllable(tf, method);
// controllable = true
</pre></blockquote>
</html>"));
    end isControllable;

    encapsulated function isObservable
      "Check oberservability of a transfer function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

        input TransferFunction tf "transfer function of a system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;
        output Boolean observable;
    protected
        StateSpace ss=StateSpace(tf);

    algorithm
        observable := StateSpace.Analysis.isObservable(ss=ss, method=method);

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
observable = TransferFunction.Analysis.<b>isObservable</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function TransferFunction.Analysis.<b>isObservable</b> checks the observability of a transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isObservable\">StateSpace.Analysis.isObservable</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2 + 2*s +1);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean controllable;

<b>algorithm</b>
  controllable := Modelica_LinearSystems2.StateSpace.Analysis.isObservable(tf, method);
// controllable = false
</pre></blockquote>
</html>"));
    end isObservable;

    encapsulated function isStabilizable
      "Check stabilizability of a transfer function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";

      output Boolean stabilizable;
    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      stabilizable := StateSpace.Analysis.isStabilizable(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
stabilizable = TransferFunction.Analysis.<b>isStabilizable</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function TransferFunction.Analysis.<b>isStabilizable</b> checks the Stabilizability of a transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isStabilizable\">StateSpace.Analysis.isStabilizable</a>.
The transfer function is stabilizable if all unstable poles are controllable.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s-1)/(s^2 - 2*s +1);

  Boolean stabilizable;

<b>algorithm</b>
   stabilizable := Modelica_LinearSystems2.TransferFunction.Analysis.isStabilizable(tf);
// stabilizable = true
</pre></blockquote>
</html>"));
    end isStabilizable;

    encapsulated function isDetectable
      "Check detectability of a transfer function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";

      output Boolean detectable;

    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      detectable := StateSpace.Analysis.isDetectable(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
detectable = TransferFunction.Analysis.<b>isDetectable</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function TransferFunction.Analysis.<b>isDetectable</b> checks the Detectability of a transfer function. Therefore, the transfer function is converted into a state space representation which is applied to <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isDetectable\">StateSpace.Analysis.isDetectable</a>.
The transfer function is detectable if all unstable poles are observable.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s-1)/(s^2 - 2*s +1);

  Boolean detectable;

<b>algorithm</b>
  detectable := Modelica_LinearSystems2.TransferFunction.Analysis.isDetectable(tf);
// detectable = false
</pre></blockquote>
</html>"));
    end isDetectable;

    encapsulated function controllabilityMatrix
      "Calculate the controllability matrix [B, A*B, ..., A^(n-1)*B] of a transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";
      output Real om[:,:];

    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      om := StateSpace.Analysis.controllabilityMatrix(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = TransferFunction.Analysis.<b>controllabilityMatrix</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the controllability matrix
<blockquote><pre>
<b>Q</b> = [<b>B</b>, <b>A</b>*<b>B</b>, ..., <b>A</b>^(n-1)*<b>B</b>]
</pre></blockquote>
<p>
of the system corresponding state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
</pre>
</blockquote>
<p>
of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

  Real Q[2,2];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.TransferFunction.Analysis.controllabilityMatrix(tf);
// Q = [0, 1, 1, -1]
</pre></blockquote>
</html>"));
    end controllabilityMatrix;

    encapsulated function observabilityMatrix
      "Calculate the observability matrix of a transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf "transfer function of a system";
      output Real om[:,:];

    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      om := StateSpace.Analysis.observabilityMatrix(ss=ss);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = TransferFunction.Analysis.<b>observabilityMatrix</b>(tf, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the observability matrix
</p>
<blockquote><pre>
<b>Q</b> = [<b>C</b>; <b>C</b>*<b>A</b>; ...; <b>C</b>*<b>A</b>^(n-1)]
</pre></blockquote>
<p>
of the system corresponding state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
</pre></blockquote>
<p>
of a transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^2+s+1);

  Real Q[2,2];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.TransferFunction.Analysis.observabilityMatrix(tf);
// Q = [1, 1, -1, 0]
</pre></blockquote>
</html>"));
    end observabilityMatrix;
  end Analysis;

  encapsulated package Design
    "Package of functions to design transfer function controllers and observers"
    import Modelica;
    extends Modelica.Icons.Package;
    encapsulated function filter
      "Generate the data record of a ZerosAndPoles transfer function from a filter description"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.ZerosAndPoles;

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
        "True, if amplitude of low pass filter at f_cut where the magnitude of the filter sagged for 3dB, otherwise unmodified filter";

      output TransferFunction filter "Filter transfer function";

    protected
      ZerosAndPoles zpFilter;
    algorithm
      zpFilter := ZerosAndPoles.Design.filter(analogFilter, filterType, order, f_cut, gain, A_ripple, normalized);
      filter := ZerosAndPoles.Conversion.toTransferFunction(zpFilter);

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
filterFunction = TransferFunction.Design.<b>filter</b>(analogFilter, filterType, order, f_cut, gain, A_ripple, normalized)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a transfer function description
of low and high pass filters. For more details see also
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Literature\">[Tietze2002]</a>, pp. 815-852.
</p>
<p>
Typical frequency responses for the four supported
low pass filter types are shown in the next figure (this figure was generated
with function <a href=\"modelica://Modelica_LinearSystems2.Examples.TransferFunction.plotBodeFilter2\">Examples.TransferFunction.plotBodeFilter2</a>):
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4Filters.png\"/> </p>
<p>
The step responses of the same low pass filters are shown in the next figure,
starting from a steady state initial filter with initial input = 0.2:
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/LowPassOrder4FiltersStepResponse.png\"/> </p>
<p>
Obviously, the frequency responses give a somewhat wrong impression of the filter
characteristics: Although Butterworth and Chebyshev filters have a significantly
steeper magnitude as the CriticalDamping and Bessel filters, the step responses
of the latter ones are much better. This means for example, that a CriticalDamping
or a Bessel filter should be selected, if a filter is mainly used to make
a non-linear inverse model realizable.
</p>
<p>
Typical frequency responses for the four supported high pass filter types are shown
in the next figure:
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4Filters.png\"/> </p>
<p>
The corresponding step responses of these high pass filters are shown in the next figure:
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/HighPassOrder4FiltersStepResponse.png\"/> </p>
<p>
All filters are available in <b>normalized</b> (default) and non-normalized form.
In the normalized form, the amplitude of the filter transfer function at the cutoff
frequency is 1/sqrt(2) (= 3 dB). Note, when comparing the filters of this function
with other software systems, the setting of \"normalized\" has to be selected
appropriately. For example, the signal processing toolbox of Matlab provides
the filters in non-normalized form and therefore a comparison makes only sense,
if normalized = <b>false</b> is set.
</p>

<h4>Example</h4>
<blockquote><pre>
  Types.AnalogFilter analogFilter=Modelica_LinearSystems2.Types.AnalogFilter.CriticalDamping;
  Integer order=2;
  Modelica.SIunits.Frequency f_cut=10;

  TransferFunction tf_filter;

algorithm
  tf_filter=Modelica_LinearSystems2.TransferFunction.Design.filter(
    order=order,
    f_cut=f_cut,
    analogFilter=analogFilter);

// tf_filter = 9530.93/(s^2 + 195.253*s + 9530.93)
</pre></blockquote>
</html>"));
    end filter;

  end Design;

  encapsulated package Plot
    "Package of functions to plot transfer function responses"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.ZerosAndPoles;
    extends Modelica.Icons.Package;

  encapsulated function polesAndZeros
      "Plot poles and/or the zeros of a transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      import Modelica_LinearSystems2.Utilities.Plot;

    input TransferFunction tf "Linear system in transfer function form";
    input Boolean poles=true "= true, to plot the poles of tf" annotation(choices(checkBox=true));
    input Boolean zeros=true "= true, to plot the zeros of tf" annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
      defaultDiagram = Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

  algorithm
    StateSpace.Plot.polesAndZeros(StateSpace(tf), poles, zeros, defaultDiagram=defaultDiagram, device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
TransferFunction.Plot.<b>polesAndZeros</b>(tf);
   or
TransferFunction.Plot.<b>polesAndZeros</b>(
  tf,
  poles=true,
  zeros=true,
  plot=true,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>());
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a pole-zero-map of the poles and zeros of a transfer function.
The Boolean inputs \"poles\" and \"zeros\" define what to plot. If Boolean input
\"plot = true\", the pole-zero-map is plotted. If false, only the diagram is generated
and returned as output argument. The records \"defaultDiagram\" and \"device\" allow
to set various layout options and the size and location of the diagram on the screen.
</p>

<h4>Example</h4>
<p>
The example <a href=\"modelica://Modelica_LinearSystems2.Examples.TransferFunction.plotPolesAndZeros\">
Modelica_LinearSystems2.Examples.TransferFunction.plotPolesAndZeros</a>
defines two transfer functions. The second one is defined as:
</p>
<blockquote><pre>
TransferFunction s   = TransferFunction.s();
TransferFunction tf2 = (s^3 + 4*s + 1)/(s^4 + 2*s^3 + 3*s^2 + 4*s);

Modelica_LinearSystems2.TransferFunction.Plot.polesAndZeros(
  tf=tf2,
  defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
    heading=\"Poles and zeros of \" + String(tf2)));
</pre></blockquote>
<p>
and results in
</p>
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/PolesAndZerosTF.png\">
</blockquote>
</html>"));
  end polesAndZeros;

    encapsulated function bode "Plot transfer function as bode plot"
      import Modelica;
      import Modelica.Utilities.Strings;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Plot;
      import SI = Modelica.SIunits;

      input TransferFunction tf "Transfer function to be plotted";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange=true
        "True, if abszissa range is automatically determined";
      input SI.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false" annotation(Dialog(enable=not autoRange));
      input SI.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                                  annotation(Dialog(enable=not autoRange));

      input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
      input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot: "
             + String(tf)));

      input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));

    protected
      SI.AngularVelocity w[nPoints];
      SI.Frequency f[nPoints];
      SI.Conversions.NonSIunits.Angle_deg phi[nPoints];
      Real A[nPoints];
      Boolean OK;
      Complex c;
      SI.Angle phi_old;
      Complex numZeros[:];
      Complex denZeros[:];

      Plot.Records.Curve curves[2];
      Integer i;
      Plot.Records.Diagram diagram2[2];

    algorithm
            // Determine frequency vector f
      if autoRange then
        (numZeros,denZeros) := TransferFunction.Analysis.zerosAndPoles(tf);
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
            denZeros);

      // Compute magnitude/phase at the frequency points
      phi_old := 0.0;
      for i in 1:nPoints loop
        w[i] := SI.Conversions.from_Hz(f[i]);
        c := TransferFunction.Analysis.evaluate(tf, Complex(0, w[i]), 1e-10);
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
TransferFunction.Plot.<b>bode</b>(tf)
   or
TransferFunction.Plot.<b>bode</b>(
  tf,
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

  encapsulated function timeResponse
      "Plot the time response of a system represented by a transfer function. The response type is selectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.TransferFunction tf;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
    input Real x0[TransferFunction.Analysis.denominatorDegree(tf)]=zeros(
      TransferFunction.Analysis.denominatorDegree(tf)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
      defaultDiagram=
      Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
      "Time response of  tf = " + String(tf)));

    protected
    Plot.Records.Curve curve;
    Plot.Records.Diagram diagram2;
    Real y[:,1,1] "Output response";
    Real t[:] "Time vector: (number of samples)";

  algorithm
    (y,t) := TransferFunction.Analysis.timeResponse(
          tf,
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
TransferFunction.Plot.<b>timeResponse</b>(tf);
   or
TransferFunction.Plot.<b>timeResponse</b>(
  tf,
  dt,
  tSpan,
  response,
  x0,
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
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);

  Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(tf, dt=0.02, tSpan=3, response=response)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/timeResponseTF.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end timeResponse;

  encapsulated function impulse "Impulse response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.TransferFunction tf "transfer function";
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

    input Real x0[TransferFunction.Analysis.denominatorDegree(tf)]=zeros(
      TransferFunction.Analysis.denominatorDegree(tf)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
      defaultDiagram=
      Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response of  tf = "
        + String(tf)));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "Type of time response";
  algorithm
    Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(
      tf=tf,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
TransferFunction.Plot.<b>impulse</b>(tf)
   or
TransferFunction.Plot.<b>impulse</b>(
  tf,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the impulse response of a transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.impulse(tf, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/impulseResponseTF.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end impulse;

  encapsulated function step "Step response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.TransferFunction tf;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";
    input Real x0[TransferFunction.Analysis.denominatorDegree(tf)]=zeros(
        TransferFunction.Analysis.denominatorDegree(tf)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  tf = "
           + String(tf)));

  algorithm
    Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(
      tf=tf,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
TransferFunction.Plot.<b>step</b>(tf)
   or
TransferFunction.Plot.<b>step</b>(
  tf,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the step response of a transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf =(s + 1)/(s^2 + 5*s + 12);

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.step(tf, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/stepResponseTF.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end step;

  encapsulated function ramp "Ramp response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.TransferFunction tf;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";
    input Real x0[TransferFunction.Analysis.denominatorDegree(tf)]=zeros(
        TransferFunction.Analysis.denominatorDegree(tf)) "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  tf = "
           + String(tf)));

  algorithm
    Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(
      tf=tf,
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
TransferFunction.Plot.<b>ramp</b>(tf)
   or
TransferFunction.Plot.<b>ramp</b>(
  tf,
  dt,
  tSpan,
  x0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the ramp response of a transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf =(2*s^2 + 7*s + 13)/(s^3 + 6*s^2 + 17*s + 12);

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.ramp(tf)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/rampResponseTF.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end ramp;

  encapsulated function initialResponse "Initial condition response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input Modelica_LinearSystems2.TransferFunction tf;
      input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";
      input Real y0 "Initial output (for initial condition plot)";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Initial response of  tf = " + String(tf) + "  with y0 = " + String(
            y0)));

    protected
      Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
          tf);
      Integer nx=TransferFunction.Analysis.denominatorDegree(tf);
      Real x0[nx]=Modelica.Math.Matrices.equalityLeastSquares(
            ss.A,
            fill(0, nx),
            ss.C,
            vector(y0)) "Initial state vector (for initial condition plot)";
  algorithm
    Modelica_LinearSystems2.TransferFunction.Plot.timeResponse(
      tf=tf,
      dt=dt,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
TransferFunction.Plot.<b>initialResponse</b>(tf)
   or
TransferFunction.Plot.<b>initialResponse</b>(
  tf,
  dt,
  tSpan,
  y0,
  columnLabels,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the initial response, i.e. the zeros input response of a transfer function. It is based on <a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf = (s + 1)/(s^2 + 5*s + 12);
  Real y0=1;

<b>algorithm</b>
  Modelica_LinearSystems2.TransferFunction.Plot.initialResponse(tf,y0=y0, dt=0.02, tSpan=3)
//  gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/TransferFunction/initialResponseTF.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Plot.ramp\">ramp</a>
</p>
</html>"));
  end initialResponse;

  end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of TransferFunction data record"
    import Modelica;
    extends Modelica.Icons.Package;
    encapsulated function toZerosAndPoles
      "Convert a TransferFunction into a ZerosAndPoles object"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.TransferFunction;
      import Complex;

      input TransferFunction tf "Transfer function of a system";
      output ZerosAndPoles zp(
        redeclare Real n1[ZerosAndPoles.Internal.numberOfRealZeros2(tf)],
        redeclare Real n2[integer((size(tf.n, 1) - 1 -
          ZerosAndPoles.Internal.numberOfRealZeros2(tf))/2),2],
        redeclare Real d1[ZerosAndPoles.Internal.numberOfRealPoles(tf)],
        redeclare Real d2[integer((size(tf.d, 1) - 1 -
          ZerosAndPoles.Internal.numberOfRealPoles(tf))/2),2]) "ZerosAndPoles object";
    protected
      Complex z[:];
      Complex p[:];
      Real k;
    algorithm
      (z,p,k) := TransferFunction.Analysis.zerosAndPoles(tf);
      zp := ZerosAndPoles(z, p, k);
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = TransferFunction.Conversion.<b>toZerosAndPoles</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a ZerosAndPoles record
</p>
<blockquote><pre>
          product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
zp = k * ---------------------------------------------------------
          product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
<p>
of a transfer function representated by numerator and denominator polynomial.
The poles and zeros and the gain <tt>k</tt> are computed
(<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">zerosAndPoles</a>)
and are used as inputs the ZerosAndPoles constructor.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction dtf = 1/(s^2 + 3*s +2)

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles(tf);
//  zp = 1/( (s + 1)*(s + 2) )
</pre></blockquote>
</html>"));
    end toZerosAndPoles;

    function toStateSpace
      "Convert a TransferFunction into a StateSpace representation"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica.Math.Vectors;

     input Modelica_LinearSystems2.TransferFunction tf
        "Transfer function of a system";
     output Modelica_LinearSystems2.StateSpace ss(
       redeclare Real A[TransferFunction.Analysis.denominatorDegree(tf),TransferFunction.Analysis.denominatorDegree(tf)],
       redeclare Real B[TransferFunction.Analysis.denominatorDegree(tf),1],
       redeclare Real C[1,TransferFunction.Analysis.denominatorDegree(tf)],
       redeclare Real D[1,1]) "Transfer function in StateSpace form";

    protected
     Integer na=TransferFunction.Analysis.denominatorDegree(tf) + 1;
     Integer nb=TransferFunction.Analysis.numeratorDegree(tf) + 1;
     Integer nx=na - 1;
     Real a[na]=Vectors.reverse(tf.d) "Reverse element order of tf.a";
     Real b[na]=vector([Vectors.reverse(tf.n); zeros(na - nb, 1)]);
     Real d=b[na]/a[na];
    algorithm
     if nx == 0 then
       ss.A := fill(0, 0, nx);
       ss.B := fill(0, 0, 1);
       ss.C := fill(0, 1, 0);
     else
       ss.A[1:nx - 1, 1:nx] := [zeros(nx - 1, 1),identity(nx - 1)];
       ss.A[nx, 1:nx] := -a[1:na - 1]/a[na];
       ss.B := [zeros(nx - 1, 1); 1/a[na]];
       ss.C := {b[1:nx] - d*a[1:nx]};

    end if;

      ss.D := [d];

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = TransferFunction.Conversion.<b>toStateSpace</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Transforms a transfer function into state space representation.
There are an infinite number of possible realizations.
Here, the transfer function is transformed into
controller canonical form, i.e. the transfer function
</p>
<blockquote><pre>
     b4*s^4 + b3*s^3 + b2*s^2 + b1*s + b0
y = -------------------------------------- * u
     a4*s^4 + a3*s^3 + a2*s^2 + a1*s + a0
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
<b>x</b> is built up of y and of all derivatives of y up to nx-1
(nx is the dimension of the state vector):
</p>
<blockquote><pre>
<b>x</b> = {y, dy/dt, d^2y/dt^2, ..., d^(n-1)y/dt^(n-1)};
</pre></blockquote>
<p>
Note, the state vector <b>x</b> of Modelica.Blocks.Continuous.TransferFunction
is defined slightly differently.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^3 + s^2 + s +1);

<b>algorithm</b>
  ss := Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace(tf);
// ss.A = [0, 1, 0; 0, 0, 1; -1, -1, -1],
// ss.B = [0; 0; 1],
// ss.C = [1, 1, 0],
// ss.D = [0],
</pre></blockquote>
</html>"));
    end toStateSpace;

    encapsulated function toZerosAndPolesMIMO
      "Convert a TransferFunction into a zeros-and-poles representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.TransferFunction;

      input TransferFunction tf[:,:] "transfer function of a system";

      output ZerosAndPoles zp[size(tf, 1),size(tf, 2)];

    protected
      Integer ny=size(tf, 1);
      Integer nu=size(tf, 2);

    algorithm
      for iy in 1:ny loop
        for iu in 1:nu loop
          zp[iy, iu] := ZerosAndPoles(tf[iy, iu]);
        end for;
      end for;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = TransferFunction.Conversion.<b>toZerosAndPolesMIMO</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Converts a matrix of transfer functions denoted as rational polynomial function into a matrix of transfer functions in zeros-and-poles representation. The function repetitively uses
<a href=\"modelica://Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles\">toZerosAndPoles</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf = [1/(s^2 + 3*s +2);s/(s^2 + 2*s +1)]

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.TransferFunction.Conversion.toZerosAndPoles(tf);
//  zp = [1/( (s + 1)*(s + 2) ); s/( (s + 1)^2 )]
</pre></blockquote>
</html>"));
    end toZerosAndPolesMIMO;

    function toMatrices
      "Convert a TransferFunction into the matrices A, B, C of a StateSpace"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica.Math.Vectors;

      input Modelica_LinearSystems2.TransferFunction tf
        "Transfer function of a system";

      output Real ABCD[size(tf.d,1),size(tf.d,1)];

    protected
     Integer na=TransferFunction.Analysis.denominatorDegree(tf) + 1;
     Integer nb=TransferFunction.Analysis.numeratorDegree(tf) + 1;
     Integer nx=na - 1;
     Real A[nx,nx];
     Real B[nx,1];
     Real C[1,nx];
     Real D[1,1];
     Real a[na]=Vectors.reverse(tf.d) "Reverse element order of tf.a";
     Real b[na]=vector([Vectors.reverse(tf.n); zeros(na - nb, 1)]);
     Real d=b[na]/a[na];
    algorithm
     if nx == 0 then
       A := fill(
           0,
           0,
           nx);
       B := fill(
           0,
           0,
           1);
       C := fill(
           0,
           1,
           0);
     else
       A[1:nx - 1, 1:nx] := [zeros(nx - 1, 1),identity(nx - 1)];
       A[nx, 1:nx] := -a[1:na - 1]/a[na];
       B := [zeros(nx - 1, 1); 1/a[na]];
       C := {b[1:nx] - d*a[1:nx]};
    end if;
      D := [d];
      ABCD := [A,B;C,D];
     annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(A, B, C, D) = TransferFunction.Conversion.toStateSpace<b>toStateSpace</b>(tf)
</pre></blockquote>

<h4>Description</h4>
<p>
Transforms a transfer function into state space representation. The outputs are the system functions A, B, C, D.
There are an infinite number of possible realizations.
Here, the transfer function is transformed into
controller canonical form, i.e. the transfer function
</p>
<blockquote><pre>
     b4*s^4 + b3*s^3 + b2*s^2 + b1*s + b0
y = -------------------------------------- * u
     a4*s^4 + a3*s^3 + a2*s^2 + a1*s + a0
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
<b>x</b> is built up of y and of all derivatives of y up to nx-1
(nx is the dimension of the state vector):
</p>
<blockquote><pre>
<b>x</b> = {y, dy/dt, d^2y/dt^2, ..., d^(n-1)y/dt^(n-1)};
</pre></blockquote>
<p>
Note, the state vector <b>x</b> of Modelica.Blocks.Continuous.TransferFunction
is defined slightly differently.
</p>

<h4>Example</h4>
<blockquote><pre>
  TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();
  Modelica_LinearSystems2.TransferFunction tf=(s+1)/(s^3 + s^2 + s +1);

<b>algorithm</b>
  (A, B, C, D) := Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace(tf);
// A = [0, 1, 0; 0, 0, 1; -1, -1, -1],
// B = [0; 0; 1],
// C = [1, 1, 0],
// D = [0],
</pre></blockquote>
</html>"));
    end toMatrices;
  end Conversion;

  encapsulated package Import
    "Package of functions to generate a TransferFunction data record from imported data"
    import Modelica;
    extends Modelica.Icons.Package;

  encapsulated function fromFile
      "Generate a TransferFunction data record by reading numenator coefficients and denominator coefficients from a file (default file name is tf.mat)"
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.Math.Polynomial;

    input String fileName="tf.mat" "Name of the transfer function data file"   annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="transfer function data file")));
    input String numName="n" "Name of the numenator of the transfer function";
    input String denName="d" "Name of the denominator of the transfer function";

    protected
    Integer numSize[2] = Streams.readMatrixSize(fileName, numName) annotation(__Dymola_allowForSize=true);
    Integer denSize[2] = Streams.readMatrixSize(fileName, denName) annotation(__Dymola_allowForSize=true);

    Real num[numSize[1],numSize[2]]=
      Streams.readRealMatrix(fileName, numName, numSize[1], numSize[2])
      "Numenator coefficients";
    Real den[denSize[1],denSize[2]]=
      Streams.readRealMatrix(fileName, denName, denSize[1], denSize[2])
      "Denominator coefficients";
    Integer ns2=numSize[2] annotation(__Dymola_allowForSize=true);
    Integer ds2=denSize[2] annotation(__Dymola_allowForSize=true);
    public
   output TransferFunction tf(n=fill(0,ns2),d=fill(0,ds2)) "transfer function";

  algorithm
    tf.n := vector(num);
    tf.d := vector(den);
    tf.uName := numName;
    tf.yName := denName;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> TransferFunction.Import.<b>fromFile</b>(fileName, numName, denName)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Reads and loads a transfer function from a mat-file <tt>fileName</tt>. The file must contain the names of the vector with the polynomial coefficients of numerator and denominator

<h4>Example</h4>
<blockquote><pre>


<b>algorithm</b>
  tf:=Modelica_LinearSystems2.TransferFunction.Import.fromFile(\"tf.mat\", \"n\", \"d\");
//  tf = (s^2 + 2*s + 3)/(4*s^2 + 5*s + 6)
</pre></blockquote>
</html>"));
  end fromFile;

  function fromModel
    "Generate a TransferFunction data record from a state space representation resulted from linearization of a model"

    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

    input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel(translate=true)));
    input Real T_linearize = 0
        "point in time of simulation to linearize the model";
    input String fileName = "dslin" "Name of the result file";

    protected
    String fileName2 = fileName + ".mat";
    Boolean OK1 = simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
    Boolean OK2 = importInitial("dsfinal.txt");
    Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize+1);
    Integer xuy[3] = ReadSystemDimension(fileName2, "ABCD");
    Integer nx = xuy[1];
    Integer nu = xuy[2];
    Integer ny = xuy[3];
    Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(fileName2, "ABCD", nx + ny, nx + nu);
    String xuyName[nx + nu + ny] = readStringMatrix(fileName2, "xuyName", nx + nu + ny);

    StateSpace result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "Model linearized at initial point";
    public
    output TransferFunction tf[:,:];

  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    result.uNames := xuyName[nx + 1:nx + nu];
    result.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    result.xNames := xuyName[1:nx];

    tf := Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunctionMIMO(result);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tf = TransferFunction.Import.<b>fromModel</b>(modelName, T_linearize, fileName)
</pre></blockquote>

<h4>Description</h4>
<p>Generate a matrix of TransferFunction data records by linearization of a model
defined by modelName. The linearization is performed at time T_linearize of
the simulation. The system is genrated by using
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Import.fromFile\">StateSpace.Import.fromFile</a>
followed by a conversion from sate space to transfer function representation.
</p>

<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
  Real T_linearize = 5;

<b>algorithm</b>
  tf = Modelica_LinearSystems2.TransferFunction.Import.fromModel(modelName, T_linearize);

//  tf = [(0.13*s^4 + 0.05558*s^3 + 1.12241*s^2 - 5.16971*s + 9.04744)/(s^6 + 0.09*s^5 + 9.13717*s^4 - 32.0637*s^3 + 58.78*s^2 + 6.3659e-014*s - 1.1703e-014);
          (0.13*s^4 + 0.05558*s^3 + 1.12241*s^2 - 5.16971*s + 9.04744)/(s^5 + 0.09*s^4 + 9.13717*s^3 - 32.0637*s^2 + 58.78*s - 2.7929e-015);
          (-0.014*s^2 + 0.31906*s - 0.8106)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.014*s^3 + 0.31906*s^2 - 0.8106*s)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.1*s^2 - 0.160918*s - 0.21842)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78);
          (-0.1*s^3 - 0.160918*s^2 - 0.21842*s)/(s^4 + 0.09*s^3 + 9.13717*s^2 - 32.0637*s + 58.78)]
</pre></blockquote>
</html>"));
  end fromModel;

  end Import;

  encapsulated package Internal
    "Package of internal material of record TransferFunction (for advanced users only)"
    extends Modelica.Icons.InternalPackage;
    import Modelica;

    encapsulated function readLength
      "Read the number n of coefficients written in a [n,1]-matrix"
      import Modelica.Utilities.Streams;

      input String fileName = "tf.mat" "Name of the transfer function data file";
      input String polyName = "n"
        "Name of the polynominal (numenator or denominator) coefficients of the transfer function"          annotation(Dialog);
      output Integer result;
    protected
      Integer polySize[2] = Streams.readMatrixSize(fileName, polyName);

    algorithm
      result := polySize[2];
    end readLength;

    encapsulated function isControllableAndObservableSISO
      "To check whether a SISO system is controllable and observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

    input TransferFunction tf;

      output Boolean controllableAndObservable;
    protected
      StateSpace ss=StateSpace(tf);

    algorithm
      controllableAndObservable := StateSpace.Internal.isControllableAndObservableSISO(ss=ss);

    end isControllableAndObservableSISO;

  end Internal;

  annotation (defaultComponentName="transferFunction", Documentation(info="<html>
<p>
This record defines the transfer function between the input signal u
and the output signal y by the coefficients of the numerator and denominator
polynomials n(s) and d(s) respectively:
</p>
<pre>        n(s)
   y = ------ * u
        d(s)
</pre>
<p>
The order of the numerator
polynomial can be larger as the order of the denominator polynomial
(in such a case, the transfer function can not be
transformed to a StateSpace system, but other operations are possible).
</p>
<p>
Example: The transfer function
</p>
<pre>             2*s+3
   y = ----------------- * u
        4*s^2 + 5*s + 6
</pre>
<p>
is transformed in the following way in a TransferFunction record:
</p>
<pre>
   <b>import</b> Modelica_LinearSystems2.TransferFunction;
   <b>import</b> Modelica.Utilities.Streams;
   TransferFunction tf(n={2,3}, d={4,5,6});
   print(\"y = \" + TransferFunction.'String'(tf) + \" * u\");
   // prints the following string to the output window:
   //   y = (2*s + 3) / (4*s^2 + 5*s + 6) * u
</pre>
</html>"));
end TransferFunction;
