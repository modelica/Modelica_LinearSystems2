within Modelica_LinearSystems2;
operator record StateSpace
  "Continuous state space description of a linear, time invariant differential equation system (data + operations)"

  Real A[:, size(A, 1)]
    annotation (Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real B[size(A, 1), :]
    annotation (Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real C[:, size(A, 1)]
    annotation (Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real D[size(C, 1), size(B, 2)]
    annotation (Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));

  String yNames[size(C, 1)]=fill("", size(C, 1)) "Names of the output signals"
    annotation (Dialog(group="Signal names"));
  String xNames[size(A, 1)]=fill("", size(A, 1)) "Names of the states"
    annotation (Dialog(group="Signal names"));
  String uNames[size(B, 2)]=fill("", size(B, 2)) "Names of the input signals"
    annotation (Dialog(group="Signal names"));

  encapsulated operator 'constructor'
    "Collection of operators to construct a StateSpace data record"
    import Modelica;
    import Modelica_LinearSystems2;

    function fromABCDMatrices
      "Generate a StateSpace data record from A, B, C and D matrices"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input Real A[:, size(A, 1)] "System matrix";
      input Real B[size(A, 1), :]=fill(
              0.0,
              size(A, 1),
              0) "Control matrix";
      input Real C[:, size(A, 1)]=fill(
              0.0,
              0,
              size(A, 1)) "Output matrix";
      input Real D[size(C, 1), size(B, 2)]=fill(
              0.0,
              0,
              0) "Feed-forward matrix";

      input String uNames[size(B, 2)]=fill("", size(B, 2)) "Names of the input signals";
      input String yNames[size(C, 1)]=fill("", size(C, 1)) "Names of the output signals";
      input String xNames[size(A, 2)]=fill("", size(A, 2)) "Names of the states";

      output StateSpace result(
        redeclare Real A[size(A, 1), size(A, 2)],
        redeclare Real B[size(B, 1), size(B, 2)],
        redeclare Real C[size(C, 1), size(C, 2)],
        redeclare Real D[size(D, 1), size(D, 2)],
        redeclare String uNames[size(B, 2)],
        redeclare String yNames[size(C, 1)],
        redeclare String xNames[size(A, 2)]) "State space record";

    algorithm
      result.A := A;
      result.B := B;
      result.C := C;
      result.D := D;
      result.uNames := uNames;
      result.yNames := yNames;
      result.xNames := xNames;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
ss = StateSpace.&apos;constructor&apos;.<b>fromABCDMatrices</b>(A, B, C, D)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss with
</p>
<blockquote><pre>
ss.A = A;
ss.B = B;
ss.C = C;
ss.D = D;
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Real A[1,1] = [1];
  Real B[1,1] = [1];
  Real C[1,1] = [1];
  Real D[1,1] = [0];

public
  StateSpace ss;

<b>algorithm</b>
  ss := 'constructor'.fromABCDMatrices(A, B, C, D);
  // ss.A = [1]
  // ss.B = [1]
  // ss.C = [1]
  // ss.D = [0]
</pre></blockquote>
</html>"));
    end fromABCDMatrices;

    function fromReal "Generate a StateSpace data record from a real value"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

      input Real r "Value of real variable";
      output StateSpace ss(
        redeclare Real A[0, 0],
        redeclare Real B[0, 1],
        redeclare Real C[1, 0],
        redeclare Real D[1, 1]) "= r";

    algorithm
      ss.D[1, 1] := r;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
ss = StateSpace.&apos;constructor&apos;.<b>fromReal</b>(r)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss from a real value, i.e. a state-space system without a state and an output without dynamics:
</p>
<blockquote><pre>
y = r*u
</pre></blockquote>
<p>
Therefore, the matrices are defined by
</p>
<blockquote><pre>
ss.A = fill(0,0,0);
ss.B = fill(0,0,1);
ss.C = fill(0,1,0);
ss.D = [r];
</pre></blockquote>
</html>"));
    end fromReal;

    function fromTransferFunction =
      Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace
      "Generate a StateSpace data record from a transfer function" annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
ss = StateSpace.&apos;constructor&apos;.<b>fromTransferFunction</b>(tf)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss from a transfer function tf.
For the simplicity of implementation, this function directly extends from
<a href=\"Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace\">TransferFunction.Conversion.toStateSpace</a>.
</p>
</html>"));

    function fromZerosAndPoles =
      Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace
      "Generate a StateSpace data record from a zeros-and-poles system"
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote>
<pre>
ss = StateSpace.&apos;constructor&apos;.<b>fromZerosAndPoles</b>(zp)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss from a zeros-poles-gain system zp.
For the simplicity of implementation, this function directly extends from
<a href=\"Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace\">ZerosAndPoles.Conversion.toStateSpace</a>.
</p>
</html>"));

    annotation (Documentation(info="<html>
<p>This package contains the default constructors for a data record of state space system. </p>
</html>"), Icon(graphics={Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25.0), Rectangle(
            lineColor={128,128,128},
            fillPattern=FillPattern.None,
            extent={{-100,-100},{100,100}},
            radius=25.0)}));
  end 'constructor';

  encapsulated operator '-'
    "Collection of operators for subtraction of state space systems"
    import Modelica;

    function subtract
      "Subtraction of two state space systems connected in parallel (= inputs are the same, outputs of the two systems are subtracted)"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss1 "State-space system 1";
      input StateSpace ss2 "State-space system 2 is subtracted from system 1";
      output StateSpace result(
        redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1), size(ss1.A, 2) + size(
          ss2.A, 2)],
        redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1), size(ss1.B, 2)],
        redeclare Real C[size(ss1.C, 1), size(ss1.C, 2) + size(ss2.C, 2)],
        redeclare Real D[size(ss1.D, 1), size(ss1.D, 2)]) "= ss1 - ss2";
    protected
      Integer nx1=size(ss1.A, 1);
      Integer nx2=size(ss2.A, 1);
    algorithm
      result.A := [ss1.A, zeros(nx1, nx2); zeros(nx2, nx1), ss2.A];
      result.B := [ss1.B; ss2.B];
      result.C := [ss1.C, -ss2.C];
      result.D := ss1.D - ss2.D;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.&apos;-&apos;.<b>subtract</b>(ss1, ss2)
</pre></blockquote>

<h4>Description</h4>
<p>
This operator function computes the subtraction of two state space systems connected in parallel,
i.e. the inputs are the same and the outputs of the two systems are subtracted.
Therefore, the systems must have the same number of inputs and outputs but not the same
number of states. The resulting system has an order of system_order1 + system_order2.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
ss3 := ss1 - ss2;
</pre> </blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss1 = StateSpace(A=[-1, 0; 0, -2], B=[1; 2], C=[0, 1], D=[0]);
  StateSpace ss2 = StateSpace(A=[-3, 0; 0, -4], B=[3; 4], C=[0, 2], D=[0]);

  StateSpace ss3;

<b>algorithm</b>
  ss3 := ss1 - ss2;
// ss.A = [-1, 0, 0, 0; 0, -2, 0, 0; 0, 0, -3, 0; 0, 0, 0, -4],
// ss.B = [1; 2; 3; 4],
// ss.C = [0, 1, 0, -2],

// ss.D = [0],
</pre></blockquote>
</html>"));
    end subtract;

    function negate
      "Unary minus (state space system where the output is multiplied by a gain of -1)"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State-space system";
      output StateSpace result(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]) "= -ss";
    algorithm
      result.A := ss.A;
      result.B := ss.B;
      result.C := -ss.C;
      result.D := -ss.D;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.&apos;-&apos;.<b>negate</b>(ss1)
</pre></blockquote>

<h4>Description</h4>
<p>
This operator function negates the state space system, i.e. the output ss
is the negation ot the state space input ss1.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
ss := -ss1;
</pre> </blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss1 = StateSpace(A=[-1, 3; 0, -2], B=[1; 2], C=[0.2, 1], D=[0.17]);

  StateSpace ss;

<b>algorithm</b>
  ss := -ss1;
// ss.A = [-1, 3; 0, -2],
// ss.B = [1; 2],
// ss.C = [-0.2, -1],
// ss.D = [-0.17],
</pre></blockquote>
</html>"));
    end negate;
    annotation (Documentation(info="<html>
<p>
This package contains operators for subtraction of state space records.
</p>
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
            points={{-50,0},{50,0}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end '-';

  encapsulated operator function '+'
    "Parallel connection of two state space systems (= inputs are the same, outputs of the two systems are added)"
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "State space system 1";
    input StateSpace ss2
      "State space system 2 is added in parallel to system 1";
    output StateSpace result(
      redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1), size(ss1.A, 2) + size(
        ss2.A, 2)],
      redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1), size(ss1.B, 2)],
      redeclare Real C[size(ss1.C, 1), size(ss1.C, 2) + size(ss2.C, 2)],
      redeclare Real D[size(ss1.D, 1), size(ss1.D, 2)]) "= ss1 + ss2";
  protected
    Integer nx1=size(ss1.A, 1);
    Integer nx2=size(ss2.A, 1);
  algorithm
    result.A := [ss1.A, zeros(nx1, nx2); zeros(nx2, nx1), ss2.A];
    result.B := [ss1.B; ss2.B];
    result.C := [ss1.C, ss2.C];
    result.D := ss1.D + ss2.D;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.<b>&apos;+&apos;</b>(ss1, ss2)
</pre></blockquote>

<h4>Description</h4>
<p>
This operator function computes the addition of two state space systems connected in parallel,
i.e. the inputs are the same and the outputs of the two systems are added.
Therefore, the systems must have the same number of inputs and outputs but not the same
number of states. The resulting system has an order of system_order1 + system_order2.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
ss3 := ss1 + ss2;
</pre> </blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss1 = StateSpace(A=[-1, 0; 0, -2], B=[1; 2], C=[0, 1], D=[0]);
  StateSpace ss2 = StateSpace(A=[-3, 0; 0, -4], B=[3; 4], C=[0, 2], D=[0.2]);

  StateSpace ss3;

<b>algorithm</b>
  ss3 := ss1 - ss2;
// ss.A = [-1, 0, 0, 0; 0, -2, 0, 0; 0, 0, -3, 0; 0, 0, 0, -4],
// ss.B = [1; 2; 3; 4],
// ss.C = [0, 1, 0, 2],

// ss.D = [0.2],
</pre></blockquote>

</html>"));
  end '+';

  encapsulated operator function '*'
    "Series connection of two state space systems"
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "State space system 1";
    input StateSpace ss2 "State space system 2";
    output StateSpace result(
      redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1), size(ss1.A, 2) + size(
        ss2.A, 2)],
      redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1), size(ss2.B, 2)],
      redeclare Real C[size(ss1.C, 1), size(ss1.C, 2) + size(ss2.C, 2)],
      redeclare Real D[size(ss1.D, 1), size(ss2.D, 2)])
      "y = G(s)*u = G(ss1)*G(ss2)*u";
  protected
    Integer nx1=size(ss1.A, 1);
    Integer nx2=size(ss2.A, 1);
  algorithm
    result.A := [ss1.A, ss1.B*ss2.C; zeros(nx2, nx1), ss2.A];
    result.B := [ss1.B*ss2.D; ss2.B];
    result.C := [ss1.C, ss1.D*ss2.C];
    result.D := ss1.D*ss2.D;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.<b>&apos;*&apos;</b>(ss1, ss2)
</pre></blockquote>

<h4>Description</h4>
<p>
This operator function computes the addition of two state space systems connected in series.
<!--,
i.e. the inputs are the same and the outputs of the two systems are added.
Therefore, the systems must have the same number of inputs and outputs but not the same
number of states. The resulting system has an order of system_order1 + system_order2.
-->
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
result := ss1 * ss2;
</pre> </blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss1 = StateSpace(A=[-1, 0; 0, -2], B=[1; 2], C=[0, 1], D=[0]);
  StateSpace ss2 = StateSpace(A=[-3, 0; 0, -4], B=[3; 4], C=[0, 2], D=[0.2]);

  StateSpace ss3;

<b>algorithm</b>
  ss3 := ss1 - ss2;
// ss.A = [-1, 0, 0, 0; 0, -2, 0, 0; 0, 0, -3, 0; 0, 0, 0, -4],
// ss.B = [0.2; 0.4; 3; 4],
// ss.C = [0, 1, 0, 0],
// ss.D = [0],
</pre></blockquote>
</html>"));
  end '*';

  encapsulated operator function '=='
    "Check whether two state space systems have identical matrices"
    import Modelica.Math.Matrices.isEqual;
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "State space system 1";
    input StateSpace ss2 "State space system 2";
    input Real eps(min=0) = 0
      "Two elements e1 and e2 of the two systems are identical if abs(e1-e2) <= eps";
    output Boolean same "=true, if the two systems are identical";
  algorithm
    same := isEqual(ss1.A, ss2.A, eps) and
      isEqual(ss1.B, ss2.B, eps) and
      isEqual(ss1.C, ss2.C, eps) and
      isEqual(ss1.D, ss2.D, eps);
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
same = StateSpace.<b>&apos;==&apos;</b>(ss1, ss2)
</pre></blockquote>

<h4>Description</h4>
<p>
This operator function returns true, if all appropriate matrices of two state space systems
ss1 and ss2 are identical. False is returned in any other case.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
same := ss1 == ss2;
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss1 = StateSpace(A=[-1, 0; 0, -2], B=[1; 2], C=[0, 1], D=[0]);
  StateSpace ss2 = StateSpace(A=[-3, 0; 0, -4], B=[3; 4], C=[0, 2], D=[0.2]);
  StateSpace ss3 = StateSpace(A=[-3, 0; 0, -4], B=[3; 4], C=[0, 2], D=[0.2]);

  ss1 == ss2;
// false

  ss2 == ss3;
// true
</pre></blockquote>
</html>"));
  end '==';

  encapsulated operator function 'String'
    "Transform state space into a String representation"
    import Modelica;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica.Utilities.Strings;

    input StateSpace ss "State space system to be transformed";
    input Integer significantDigits=12
      "Number of significant digits that are shown";
    input String name="ss" "Independent variable name used for printing";
    output String s="";

  protected
    String space=Strings.repeat(5);
    String space2=Strings.repeat(3);
    String space3=Strings.repeat(8);
    Integer nx=size(ss.A, 1);
    Integer nu=size(ss.B, 2);
    Integer ny=size(ss.C, 1);
    Integer sizeD=size(ss.D, 2);
    Integer stringMaxLength;
    Boolean xNamesExist=false;
    Boolean uNamesExist=false;
    Boolean yNamesExist=false;

  algorithm
    // If system is too large, do not print the matrices
    if size(ss.A,1) > 50 or size(ss.B, 2) > 50 or size(ss.C, 1) > 50 then
      s := "System not printed since too large (only dimensions):\n" +
           "   " + name + ".A[" + String(nx) + "," + String(nx) + "]\n" +
           "   " + name + ".B[" + String(nx) + "," + String(nu) + "]\n" +
           "   " + name + ".C[" + String(ny) + "," + String(nx) + "]\n" +
           "   " + name + ".D[" + String(ny) + "," + String(nu) + "]";
      return;
    end if;

    //Checking if name arrays are empty
    for i in 1:nx loop
      xNamesExist := xNamesExist or (ss.xNames[i] <> "");
    end for;

    for i in 1:ny loop
      yNamesExist := yNamesExist or (ss.yNames[i] <> "");
    end for;

    for i in 1:nu loop
      uNamesExist := uNamesExist or (ss.uNames[i] <> "");
    end for;
    /*
    if xNamesExist then
      Modelica.Utilities.Streams.print("xNamesExist == true");
    else
      Modelica.Utilities.Streams.print("xNamesExist == false");
    end if;
*/
    stringMaxLength := max(size(ss.xNames, 1), min(size(ss.yNames, 1), 11));

    if nx == 0 and sizeD == 0 then
      s := name + ".A = []\n  " + name + ".B = []\n   " + name +
        ".C = [] \n   " + name + ".D = []";
    else
      s := "\n" + name + ".A = \n";

      //Horizontal
      // Two alternatives when printing state names
      if xNamesExist == false then
        s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
          "x1 ";
      else
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(ss.xNames[
          1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[1]), 11)) +
          " " + Strings.substring(
            ss.xNames[1],
            1,
            min(Strings.length(ss.xNames[1]), 11));
      end if;

      for i in 2:nx loop

        //Two alternatives when printing state names

        if xNamesExist == false then
          s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
             + String(i - 1))) + "x" + String(i) + " ";
        else
          s := s + " " + Strings.repeat(significantDigits + 11 - min(
            Strings.length(ss.xNames[i - 1]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11));

        end if;

        //s := s + Strings.repeat(6) + "x" + String(i);
      end for;
      s := s + "\n";

      for i in 1:nx loop
        //Vertical
        //Two alternatives when printing state names
        if xNamesExist == false then
          s := s + space + "x" + String(i) + " ";
        else
          s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
            ss.xNames[i]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11)) + " ";
        end if;

        for j in 1:nx loop
          if ss.A[i, j] >= 0 then
            s := s + " ";
          end if;
          s := s + String(ss.A[i, j], significantDigits=significantDigits) +
            Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
            ss.A[i, j]), significantDigits=significantDigits)));
        end for;
        s := s + "\n";
      end for;
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      s := s + "\n" + name + ".B = \n";
      //Horizontal
      // Two alternatives when printing state names
      if uNamesExist == false then
        s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
          "u1 ";
      else
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(ss.uNames[
          1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[1]), 11)) +
          " " + Strings.substring(
            ss.uNames[1],
            1,
            min(Strings.length(ss.uNames[1]), 11));
      end if;

      for i in 2:nu loop
        //Two alternatives when printing state names
        if uNamesExist == false then
          s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
             + String(i - 1))) + "u" + String(i) + " ";
        else
          s := s + " " + Strings.repeat(significantDigits + 11 - min(
            Strings.length(ss.uNames[i - 1]), 11)) + Strings.substring(
              ss.uNames[i],
              1,
              min(Strings.length(ss.uNames[i]), 11));
        end if;
      end for;
      s := s + "\n";
      //s := s + Strings.repeat(6) + "x" + String(i);
      for i in 1:nx loop

        //Vertical
        //Two alternatives when printing state names
        if xNamesExist == false then
          s := s + space + "x" + String(i) + " ";
        else

          s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
            ss.xNames[i]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11)) + " ";
        end if;

        for j in 1:nu loop
          if ss.B[i, j] >= 0 then
            s := s + " ";
          end if;
          s := s + String(ss.B[i, j], significantDigits=significantDigits) +
            Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
            ss.B[i, j]), significantDigits=significantDigits)));
        end for;
        s := s + "\n";
      end for;

      //--------------------------------------------------------------------------------------------------------------------------------------------------
      s := s + "\n" + name + ".C = \n";
      //Horizontal
      // Two alternatives when printing state names
      if xNamesExist == false then
        s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
          "x1 ";
      else
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(ss.xNames[
          1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[1]), 11)) +
          " " + Strings.substring(
            ss.xNames[1],
            1,
            min(Strings.length(ss.xNames[1]), 11));
      end if;

      for i in 2:nx loop
        //Two alternatives when printing state names
        if xNamesExist == false then
          s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
             + String(i - 1))) + "x" + String(i) + " ";
        else
          s := s + " " + Strings.repeat(significantDigits + 11 - min(
            Strings.length(ss.xNames[i - 1]), 11)) + Strings.substring(
              ss.xNames[i],
              1,
              min(Strings.length(ss.xNames[i]), 11));
        end if;
      end for;
      s := s + "\n";
      //s := s + Strings.repeat(6) + "x" + String(i);

      for i in 1:ny loop
        //Vertical
        //Two alternatives when printing state names
        if yNamesExist == false then
          s := s + space + "y" + String(i) + " ";
        else
          s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
            ss.yNames[i]), 11)) + Strings.substring(
              ss.yNames[i],
              1,
              min(Strings.length(ss.yNames[i]), 11)) + " ";

        end if;

        for j in 1:nx loop
          if ss.C[i, j] >= 0 then
            s := s + " ";
          end if;
          s := s + String(ss.C[i, j], significantDigits=significantDigits) +
            Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
            ss.C[i, j]), significantDigits=significantDigits)));
        end for;
        s := s + "\n";
      end for;
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      s := s + "\n" + name + ".D = \n";
      //Horizontal
      // Two alternatives when printing state names
      if uNamesExist == false then
        s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
          "u1 ";
      else
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(ss.uNames[
          1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[1]), 11)) +
          " " + Strings.substring(
            ss.uNames[1],
            1,
            min(Strings.length(ss.uNames[1]), 11));
      end if;

      for i in 2:nu loop
        //Two alternatives when printing state names
        if uNamesExist == false then
          s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
             + String(i - 1))) + "u" + String(i) + " ";
        else
          s := s + " " + Strings.repeat(significantDigits + 11 - min(
            Strings.length(ss.uNames[i - 1]), 11)) + Strings.substring(
              ss.uNames[i],
              1,
              min(Strings.length(ss.uNames[i]), 11));
        end if;
      end for;
      s := s + "\n";
      for i in 1:ny loop
        //Vertical
        //Two alternatives when printing state names
        if yNamesExist == false then
          s := s + space + "y" + String(i) + " ";
        else
          s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
            ss.yNames[i]), 11)) + Strings.substring(
              ss.yNames[i],
              1,
              min(Strings.length(ss.yNames[i]), 11)) + " ";
        end if;

        for j in 1:nu loop
          if ss.D[i, j] >= 0 then
            s := s + " ";
          end if;
          s := s + String(ss.D[i, j], significantDigits=significantDigits) +
            Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
            ss.D[i, j]), significantDigits=significantDigits)));
        end for;
        s := s + "\n";
      end for;

    end if;

    annotation (Documentation(info="<html>
<p>
Returns a pretty formated string representing the input state space ss.
The operator is used by writing just the following command:
</p>
<blockquote><pre>
ss;
</pre></blockquote>
</html>"));
  end 'String';

  encapsulated package Analysis
    "Collection of functions to analyse state space system represented by a StateSpace record"
    import Modelica;
    extends Modelica.Icons.Package;

    function analysis
      "Perform a system analysis based on the poles and zeros of the system"

      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Eigenvalue;
      import Complex;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.Utilities.Plot;
      import DymolaCommands;

      input StateSpace ss;

      input Internal.AnalyseOptions analyseOptions=
          Modelica_LinearSystems2.Internal.AnalyseOptions(
              plotEigenValues=true,
              plotInvariantZeros=true,
              plotStepResponse=true,
              plotFrequencyResponse=true,
              printSystem=true,
              printEigenValues=true,
              printEigenValueProperties=true,
              printInvariantZeros=true,
              printControllability=true,
              printObservability=true,
              headingEigenValues="Eigenvalues",
              headingInvariantzeros="Invariant zeros",
              headingStepResponse="Step response",
              headingFrequencyResponse="Frequency response",
              dB_w = false);
      input String fileName="systemReport.html"
        "Name of html-file that contains eigenvalue table";
      input String systemName=""
        "Name of system (used as heading in html file)";
      input String description="" "Description of system (used in html file)";
    protected
      String dummyFileName="dummy" + fileName;
    public
      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

    protected
      StateSpace ssBalanced = StateSpace.Transformation.toBalancedForm(ss);
      Complex j = Modelica.ComplexMath.j;
      Eigenvalue ev[size(ss.A, 1)];
      Integer nx=size(ss.A, 1);
      Integer window=0;
      Real eval[nx, 2];
      Real revec[nx, nx];
      Real levec[nx, nx];
      Complex cev[size(ss.A, 1)];
      Complex systemZeros[:]=
          Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ssBalanced);
      Boolean isStable;
      Boolean isControllable;
      Boolean isStabilizable;
      Boolean isObservable;
      Boolean isDetectable;

      Real abs_evec[nx];
      String xNames2[nx];
      String heading="Eigenvalues" "Eigen values of system";
      Eigenvalue evSorted[size(ss.A, 1)];
      Integer evIndex[size(ss.A, 1)];
      Complex zerosSorted[:];
      Integer zerosIndex[size(systemZeros, 1)];
      Integer nReal;

      Integer i;
      Integer k;
      Complex evecComplex[size(ss.A, 1), size(ss.A, 1)];
      Plot.Records.Curve curves[2];
      Plot.Records.Diagram diagram2;
      Boolean instableZeros=false;

      String filePathOnly "NOT USED: Path to fileName";
      String fileNameOnly "Name of fileName without extension and path";
      String fileExtOnly "Extension of fileName";
      String fileNameImg "General name (without extension) of file for a plot";
      String fileNameImg2="none" "Current name of file for a plot";

      Internal.AnalyseOptions analyseOptions2=analyseOptions;
    algorithm
      // ---------------------------------------------------------------------------------------------------
      // The correct HTML format generated with this function can be checked with following commands:
      //   Modelica_LinearSystems2.StateSpace.Analysis.analysis(Modelica_LinearSystems2.StateSpace(A=[2,1,1;1,1,1;1,2,2], B=[1;2.2;3], C=[2,4,6;3,8,5], D=[6;4], yNames={"y1_test","y2_test"}, xNames={"xx1","xx2","xx3"}, uNames={"u1_test"}), description="Test file in HTML format from function 'analysis'.");
      //   Modelica_LinearSystems2.StateSpace.Analysis.analysis(Modelica_LinearSystems2.StateSpace(A=[2, 1.43, 12, 3; 1, 1, 1, 43; 1, 3, 2, 2; 1, 1, 4.2, 1.2], B=[1, 2; 2.2, 3; 3, 1; 4, 0], C=[25, 1.4, 6.3, 1; 0.3, 8, 5, 1; 1, 3, 2, 2], D=[6, 4; 4, 2; 6, 5], yNames={"y1_test","y2_te","y3_"}, xNames={"xx1","x2","xxx3","xx4"}, uNames={"u1_test","u2_test"}));
      // ---------------------------------------------------------------------------------------------------

      (filePathOnly,fileNameOnly,fileExtOnly) :=
        Modelica.Utilities.Files.splitPathName(fileName);
      fileNameImg := fileNameOnly;

      // If system has no inputs and outputs, modify analyze options that do not make sense
      if size(ss.B, 2) == 0 or size(ss.C, 1) == 0 then
        analyseOptions2.plotStepResponse := false;
        analyseOptions2.plotFrequencyResponse := false;
        analyseOptions2.printControllability := false;
        analyseOptions2.printObservability := false;
      end if;

      // If system has no states, modify analyze options that do not make sense
      if nx < 1 then
         analyseOptions2.plotEigenValues          :=false;
         analyseOptions2.plotInvariantZeros       :=false;
         analyseOptions2.printEigenValues         :=false;
         analyseOptions2.printEigenValueProperties:=false;
         analyseOptions2.printInvariantZeros      :=false;
      end if;

      // If system is too large, do not print A,B,C,D matrices
      if nx > 50 or size(ss.B, 2) > 50 or size(ss.C, 1) > 50 then
         analyseOptions2.printSystem:=false;
      end if;

      // Get eigenvalues
      // ---------------
      (eval,levec,revec) := Modelica_LinearSystems2.Math.Matrices.eigenValues(
        ss.A);

      for i in 1:nx loop
        cev[i].re := eval[i, 1];
        cev[i].im := eval[i, 2];
        ev[i].ev := cev[i];
      end for;

      (evSorted,evIndex) := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

      // Build x names
      // -------------
      if size(ss.xNames, 1) <> nx or nx > 1 and ss.xNames[1]=="" then
        for i in 1:nx loop
          xNames2[i] := "x[" + String(i) + "]";
        end for;
      else
        xNames2 := ss.xNames;
      end if;

      // Whole system checks
      // ===================
      // Stability check
      isStable := true;
      for i in 1:nx loop
        isStable := isStable and ev[i].ev.re < 0;
      end for;

      if size(ss.B, 2) > 0 and size(ss.C, 1) > 0 then
        // Controllability check, stabilizability check
        isControllable := StateSpace.Analysis.isControllable(ssBalanced);
        isStabilizable := StateSpace.Analysis.isStabilizable(ssBalanced);

        // Observability check, detectability check
        isObservable := StateSpace.Analysis.isObservable(ssBalanced);
        isDetectable := StateSpace.Analysis.isDetectable(ssBalanced);
      else
        isControllable := false;
        isStabilizable := false;
        isObservable := false;
        isDetectable := false;
      end if;

      // Analysis of single eingenvalues
      ev := StateSpace.Internal.characterizeEigenvalue(ss, ev);

      // Sort eigen values according to smallest imaginary value and restore the original order
      evSorted := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

      // Analysis file
      // -------------
      Modelica.Utilities.Files.removeFile(fileName);
      Modelica.Utilities.Files.removeFile(dummyFileName);

      // Text should be printed into new file in HTML environment
      // --------------------------------------------------------
      StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
      StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);

      if analyseOptions2.printSystem then
        printSystem(
              ss,
              fileName,
              systemName,
              description);

        printSystem(
              ss,
              dummyFileName,
              systemName,
              description);
      end if;

      printHead1(
            ss,
            isStable,
            isControllable,
            isStabilizable,
            isObservable,
            isDetectable,
            fileName,
            analyseOptions=analyseOptions2);
      printHead1(
            ss,
            isStable,
            isControllable,
            isStabilizable,
            isObservable,
            isDetectable,
            dummyFileName,
            analyseOptions=analyseOptions2);
      StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);

      Modelica.Utilities.Streams.readFile(dummyFileName);

      // Plot step response
      // ------------------
      if analyseOptions2.plotStepResponse then
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print(
          "<html>\n<body>\n<p>\n<b>Step responses</b>\n</p>\n</body>\n</html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);
        StateSpace.Plot.step(ss=ssBalanced);
        fileNameImg2 := fileNameImg + "Step.png";
        DymolaCommands.Plot.ExportPlotAsImage(fileName=fileNameImg2, id=-1, includeInLog=false);
        print("<p>\n<img src=\"" + fileNameImg2 + "\">\n</p>", fileName);
      end if;

      // Plot Bode plots
      if analyseOptions2.plotFrequencyResponse then
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print("<html>\n<body>\n<p>\n<b>Bode plots</b>\n</p>\n</body>\n</html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);
        StateSpace.Plot.bodeMIMO(ss=ss, Hz=not analyseOptions.dB_w, dB=analyseOptions.dB_w);
        //     fileNameImg2 := fileNameImg + "BodeMIMO1.png";
        //     ExportPlotAsImage(
        //       fileName=fileNameImg2,
        //       id=-1,
        //       includeInLog=false);
      end if;

      // Calculate the number of real eigenvalues
      nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(cev);

      // Construct complex eigenvector matrix
      i := 1;
      while i <= nx loop
        if eval[i, 2] == 0.0 then

          for jj in 1:nx loop
            evecComplex[jj, i] := revec[jj, i] + 0*j;
          end for;
          i := i + 1;
        else
          for jj in 1:nx loop
            evecComplex[jj, i] := revec[jj, i] + revec[jj, i + 1]*j;
            evecComplex[jj, i + 1] := revec[jj, i] - revec[jj, i + 1]*j;
          end for;
          i := i + 2;
        end if;
      end while;

      if analyseOptions2.printEigenValues then
        printHead2a(
              fileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal > 0));
        if nReal > 0 then
          printTab1(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                fileName,
                analyseOptions=analyseOptions2);
        end if;

        Modelica.Utilities.Files.removeFile(dummyFileName);
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead2a(
              dummyFileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal > 0));
        if nReal > 0 then
          printTab1(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                dummyFileName,
                analyseOptions=analyseOptions2);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
        Modelica.Utilities.Streams.readFile(dummyFileName);

        printHead2b(
              fileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal < nx));
        if nReal < nx then
          printTab2(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                fileName,
                analyseOptions=analyseOptions2);
        end if;

        Modelica.Utilities.Files.removeFile(dummyFileName);
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead2b(
              dummyFileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal < nx));
        if nReal < nx then
          printTab2(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                dummyFileName,
                analyseOptions=analyseOptions2);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
        Modelica.Utilities.Streams.readFile(dummyFileName);

       // Plot eigenvalues and invariant zeros
       if analyseOptions2.plotEigenValues or
          analyseOptions2.plotInvariantZeros and size(systemZeros, 1) > 0 then
          i := 0;
          if analyseOptions2.plotEigenValues then
            i := i + 1;
            curves[i] := Plot.Records.Curve(
                  x=eval[:, 1],
                  y=eval[:, 2],
                  legend="poles",
                  autoLine=false,
                  linePattern=Plot.Types.LinePattern.None,
                  lineSymbol=Plot.Types.PointSymbol.Cross);
          end if;

          // Plot invariant zeros
          if size(systemZeros, 1) > 0 and analyseOptions2.plotInvariantZeros then
            i := i + 1;
            curves[i] := Plot.Records.Curve(
                  x=systemZeros[:].re,
                  y=systemZeros[:].im,
                  legend="zeros",
                  autoLine=false,
                  linePattern=Plot.Types.LinePattern.None,
                  lineSymbol=Plot.Types.PointSymbol.Circle);
          end if;

          diagram2 := defaultDiagram;
          diagram2.curve := curves[1:i];
          Plot.diagram(diagram2, device);
        end if;

        if analyseOptions2.printEigenValueProperties then
          printHead3(fileName);
          printTab3(
                evSorted,
                evecComplex,
                evIndex,
                cev,
                nReal,
                xNames2,
                fileName);

          Modelica.Utilities.Files.removeFile(dummyFileName);
          StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
          printHead3(dummyFileName);
          printTab3(
                evSorted,
                evecComplex,
                evIndex,
                cev,
                nReal,
                xNames2,
                dummyFileName);
          StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
          Modelica.Utilities.Streams.readFile(dummyFileName);

        end if;
      end if;

      // ZEROS
      (zerosSorted,zerosIndex) :=Modelica.ComplexMath.Vectors.sort(systemZeros);
      nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(zerosSorted);

      if analyseOptions2.printInvariantZeros then
        printHead4(fileName, printTable=(size(systemZeros, 1) > 0));
        if size(systemZeros, 1) > 0 then
          Modelica_LinearSystems2.StateSpace.Analysis.analysis.printTab4(
                zerosSorted,
                zerosIndex,
                nReal,
                fileName);
        end if;

        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead4(dummyFileName, printTable=(size(systemZeros, 1) > 0));
        if size(systemZeros, 1) > 0 then
          printTab4(
                zerosSorted,
                zerosIndex,
                nReal,
                dummyFileName);
        end if;
        k := 0;
        for i in 1:size(systemZeros, 1) loop
          if systemZeros[i].re > 0 then
            k := k + 1;
          end if;
        end for;
        if k > 0 then
          print("<p>\n<b>Note, that the system has " + String(k) +
            " zeros in the right complex half-plane.</b>\n</p>", fileName);
          print("<p>\n<b>Note, that the system has " + String(k) +
            " zeros in the right complex half-plane.</b>\n</p>", dummyFileName);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);

      end if;
      Modelica.Utilities.Streams.readFile(dummyFileName);
      Modelica.Utilities.Files.removeFile(dummyFileName);

      print("\n\nAnalysis results have been written to file \"" +
        Modelica.Utilities.Files.fullPathName(fileName) + "\"");

      // Last print of HTML environment
      // --------------------------------------------------------
      StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);

      // SUB FUNCTIONS

    public
      encapsulated function printSystem
        "Print the state space system in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;
        import Modelica_LinearSystems2;

        input StateSpace ss "State space system to analyze";
        input String fileName="systemAnalysis.html"
          "File on which the state space system is written in html format";
        input String systemName="State Space System"
          "Name of the state space system";
        input String description="" "Description of system (used in html file)";
        input String format=".6g" "Format of numbers (e.g. \"20.8e\")";

        input Boolean htmlBasics=false
          "True, if text should be printed within 'html' and 'body' environment, otherwise text printed into existing file fileName"
          annotation (Dialog(group="HTML format"));
        input Integer hSize(
          min=1,
          max=5) = 1
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));
      protected
        Integer nx=size(ss.A, 1);
        Integer nu=size(ss.B, 2);
        Integer ny=size(ss.C, 1);
        Integer c1=integer(ceil(nx/2) - 1);
        Integer c2=integer(ceil(ny/2) - 1);
        Integer dist=2;
        Boolean centered=true
          "True, if matrices columns should be centered, otherwise right aligned";

        String td_align=if centered then "  <td style=\"text-align:center\">"
             else "  <td style=\"text-align:right\">";

      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 4 then 4 else hSize;
        String heading="h" + String(hSizeOK);
        String heading2="h" + String(hSizeOK + 1);
        String heading3="h" + String(hSizeOK + 2);
        Boolean printIndices;

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(Modelica_LinearSystems2.StateSpace(A=[2,1,1;1,1,1;1,2,2], B=[1;2.2;3], C=[2,4,6;3,8,5], D=[6;4], yNames={"y1_test","y2_test"}, xNames={"xx1","xx2","xx3"}, uNames={"u1_test"}));
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(Modelica_LinearSystems2.StateSpace(A=[2, 1.43, 12, 3; 1, 1, 1, 43; 1, 3, 2, 2; 1, 1, 4.2, 1.2], B=[1, 2; 2.2, 3; 3, 1; 4, 0], C=[25, 1.4, 6.3, 1; 0.3, 8, 5, 1; 1, 3, 2, 2], D=[6, 4; 4, 2; 6, 5], yNames={"y1_test","y2_te","y3_"}, xNames={"xx1","x2","xxx3","xx4"}, uNames={"u1_test","u2_test"}));
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(ss=Modelica_LinearSystems2.StateSpace.Import.fromModel("Modelica.Mechanics.Rotational.Examples.First"), description="Test file in HTML format from function printSystem.");
        // ---------------------------------------------------------------------------------------------------

        if htmlBasics then
          // Text should be printed into new file in HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
        end if;

        print("<" + heading + ">System report</" + heading + ">", fileName);
        print("\n<" + heading2 + ">General information</" + heading2 + ">",
          fileName);
        //print("<h1>System report</h1>", fileName);
        //print("\n<h2>General information</h2>", fileName);

        if systemName == "" then
        else
          print("\n<" + heading3 + ">System name</" + heading3 + ">", fileName);
          //print("\n<h3>System name</h3>", fileName);
          print("<p>\n" + systemName + "\n</p>", fileName);
        end if;

        if description == "" then
        else
          print("\n<" + heading3 + ">Description</" + heading3 + ">", fileName);
          //print("\n<h3>Description</h3>", fileName);
          print("<p>\n" + description + "\n</p>", fileName);
        end if;

        print("\n<" + heading3 + ">Matrices</" + heading3 + ">", fileName);
        //print("\n<h3>Matrices</h3>", fileName);
        print(
          "<p>\nThe system described in the state space representation\n</p>",
          fileName);
        print(
          "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin: 20px 0 20px 20px\" "
           + "cellpadding=\"3\" border=\"0\"> ", fileName);
        print("<tr><td>der(x) </td> <td>=</td> <td> Ax</td> <td> +</td><td> Bu</td></tr>
         <tr><td> y </td>     <td>=</td> <td> Cx</td> <td> + </td><td>Du</td></tr>",
          fileName);
        print("</table>\n<p>\nis defined by\n</p>", fileName);

        // ===============================
        // Print signal names and matrices (print row and column indices if at least one matrix has more as 5 elements)
        // ===============================
        printIndices := size(ss.A, 1) > 5 or size(ss.B, 2) > 5 or size(ss.C, 1)
           > 5;
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.uNames,
                "uNames",
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.yNames,
                "yNames",
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.xNames,
                "xNames",
                fileName=fileName,
                printIndices=printIndices);

        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.A,
                "A",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.B,
                "B",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.C,
                "C",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.D,
                "D",
                format=format,
                fileName=fileName,
                printIndices=printIndices);

        if ny == 0 and nu == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has neither inputs nor outputs (and therefore matrices B, C, and D are empty matrices)!\n</p>",
            fileName);
        elseif ny == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has no outputs (and therefore matrices C and D are empty matrices)!\n</p>",
            fileName);
        elseif nu == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has no inputs (and therefore matrices B and D are empty matrices)!\n</p>",
            fileName);
        end if;

        if htmlBasics then
          // Last print of HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);
        end if;

      end printSystem;

      encapsulated function printHead1
        "Print the heading of document for characteristics in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input StateSpace ss;
        // This could be deleted sinc not used. But for reasons of beackward compatibility it is still here.
        input Boolean isStable;
        input Boolean isControllable;
        input Boolean isStabilizable;
        input Boolean isObservable;
        input Boolean isDetectable;
        input String fileName="systemHead1.html"
          "File on which the information is written in html format";

        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

        input Boolean htmlBasics=false
          "True, if text should be printed within 'html' and 'body' environment, otherwise text printed into existing file fileName"
          annotation (Dialog(group="HTML format"));
        input Integer hSize(
          min=1,
          max=5) = 2
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));

      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 5 then 5 else hSize;
        String heading="h" + String(hSizeOK);

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead1(Modelica_LinearSystems2.StateSpace(A=[2], B=[1], C=[1], D=[1]), false, false, false, true, false, htmlBasics=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        if htmlBasics then
          // Text should be printed into new file in HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
        end if;

        print("\n<" + heading + ">Characteristics</" + heading +
          ">\n<p>\nThe system\n</p>\n<p> is ", fileName);

        if analyseOptions.printControllability and analyseOptions.printObservability then
          print((if isStable then " " else "<b>not</b> ") + "stable" + "\n<br>"
             + (if isStable then if isControllable then "and it is " else
            "but it is <b>not</b> " else if isControllable then "but it is "
             else "and it is <b>not</b> ") + "controllable" + (if isStable
             then "" else "\n<br>" + (if isControllable then
            " and therefore it is " else if isStabilizable then " but it is "
             else "and is <b>not</b> ") + "stabilizable.") +
            "\n<br> The system is " + (if isObservable then " " else
            "<b>not</b> ") + "observable" + (if isStable then "" else "\n<br>"
             + (if isObservable then " and therefore it is " else if
            isDetectable then " but it is " else "and is <b>not</b> ") +
            "detectable.") + "\n<br>", fileName);
        elseif not analyseOptions.printObservability and analyseOptions.printControllability then
          print((if isStable then " " else "<b>not</b> ") + "stable" + "\n<br>"
             + (if isStable then if isControllable then "and it is " else
            "but it is <b>not</b> " else if isControllable then "but it is "
             else "and it is <b>not</b> ") + "controllable" + (if isStable
             then "" else "\n<br>" + (if isControllable then
            " and therefore it is " else if isStabilizable then " but it is "
             else "and is <b>not</b> ") + "stabilizable.") + "\n<br>", fileName);
        elseif not analyseOptions.printControllability and analyseOptions.printObservability then
          print((if isStable then " " else "<b>not</b> ") + "stable." +
            "\n<br> The system is " + (if isObservable then " " else
            "<b>not</b> ") + "observable" + (if isStable then "" else "\n<br>"
             + (if isObservable then " and therefore it is " else if
            isDetectable then " but it is " else "and is <b>not</b> ") +
            "detectable.") + "\n<br>", fileName);
        else
          print((if isStable then " " else "<b>not</b> ") + "stable." +
            "\n<br>", fileName);
        end if;

        print("</p>", fileName);

        if htmlBasics then
          // Last print of HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);
        end if;

      end printHead1;

      encapsulated function printHead2a
        "Print the heading of document for eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input String fileName="systemHead2a.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

        input Boolean printTable=true
          "True, if the system has real eigenvalues to be printed in table";
        input Integer hSize(
          min=1,
          max=5) = 3
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));
      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 5 then 5 else hSize;
        String heading="h" + String(hSizeOK);

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead2a(htmlBasics=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        print("\n<" + heading + ">Eigenvalues analysis</" + heading + ">",
          fileName);
        //print("<p>\n<b>Real eigenvalues</b>\n</p>", fileName);

        if printTable then
          print("<p>The system has the following real eigenvalues.</p>",
            fileName);
          print(
            "<table style=\"background-color:rgb(100, 100, 100); margin: 20px 0 20px 20px\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Real eigenvalues</caption>", fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> eigenvalue </td>\n  <td> T [s] </td>\n  <td> characteristics </td>",
            fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td> contribution to states</td>", fileName);
          end if;

          print("</tr>", fileName);
        else
          print("<p>\nThe system has no real eigenvalues.\n</p>", fileName);
        end if;

      end printHead2a;

      encapsulated function printHead2b
        "Print the heading of document for conjugated complex pairs in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

        input String fileName="systemHead2b.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");
        input Boolean printTable=true
          "True, if the system has complex pairs to be printed in table";

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead2b();
        // ---------------------------------------------------------------------------------------------------

        if printTable then
          print(
            "<p>The system has the following complex conjugate pairs of eigenvalues.<br>&nbsp;</p>",
            fileName);
          print(
            "<table style=\"background-color:rgb(100, 100, 100); margin: 20px 0 20px 20px\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Complex conjugate pairs of eigenvalues</caption>",
            fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> eigenvalue </td>\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n  <td> characteristics </td>",
            fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td> contribution to states</td>", fileName);
          end if;

          print("</tr>", fileName);
        else
          print(
            "<p>\nThe system has no complex conjugate eigenvalue pairs.\n</p>",
            fileName);
        end if;

      end printHead2b;

      encapsulated function printHead3
        "Print the heading of document for description in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

        input String fileName="systemHead3.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      algorithm
        print(
          "<p>\nIn the tables above, the column <b>contribution to states</b> lists for each eigenvalue the states to which the"
           +
          " corresponding modal state z[i] contributes most. This information is based on the"
           +
          " two largest absolute values of the corresponding right eigenvector (if the second large value"
           +
          " is less than 5&nbsp;% of the largest contribution, it is not shown). Note"
           +
          " the <b>right eigenvector</b> v<sub>j</sub> and the <b>left eigenvector</b> u<sub>j</sub> of A satisfy the"
           +
          " following relationships with regards to <b>eigenvalue</b> &lambda;<sub>j</sub>,"
           +
          " state vector x and modal state vector z (u<sub>j</sub><sup>H</sup> denotes the conjugate transpose of u<sub>j</sub>):"
           + " </p>" +
          " <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\">" +
          " <tr><td width=\"50\"></td>" +
          "\n    <td>A * v<sub>j</sub> = &lambda;<sub>j</sub> * v<sub>j</sub>; &nbsp;&nbsp;&nbsp;&nbsp;"
           +
          "         u<sub>j</sub><sup>H</sup> * A = &lambda;<sub>j</sub> * u<sub>j</sub><sup>H</sup>; &nbsp;&nbsp;&nbsp;&nbsp;"
           +
          "               x = V * z; &nbsp;&nbsp;&nbsp;&nbsp; V = [v<sub>1</sub>, v<sub>2</sub>, ...]</td>"
           + "           </tr>" + "\n</table>" + "\n<p>" +
          "\nIn the next table, for each state in the column <b>correlation to modal states</b>, the modal"
           +
          " states z[i] which contribute most to the corresponding state are summarized, that is"
           + " the state is mostly composed of these modal states." +
          "\nThis information is based on the two largest absolute values of row i of the"
           +
          " eigenvector matrix that is associated with eigenvalue i (if the second large value"
           +
          " is less than 5&nbsp;% of the largest contribution, it is not shown). This only holds"
           +
          " if the modal states z[i] are in the same order of magnitude. Otherwise, the listed modal states"
           + " might be not the most relevant ones.</p>", fileName);

        print(
          "<table style=\"background-color:rgb(100, 100, 100); margin: 20px 0 20px 20px\" "
           + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">\n" +
          "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           +
          "\n  <td> state </td>\n  <td> correlation to modal states </td>\n  <td> eigenvalue # </td>"
           +
          "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n  <td> T [s] </td>\n</tr>",
          fileName);

      end printHead3;

      encapsulated function printHead4
        "Print the heading of document for invariant zeros in html format on file"
        import Modelica;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input String fileName="systemHead4.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");
        input Boolean printTable=true
          "True, if the system has complex pairs to be printed in table";

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead4(htmlEnv=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        if printTable then
          print("<p>The system has the following invariant zeros.<br>&nbsp;</p>",
            fileName);
          print(
            "\n<table style=\"background-color:rgb(100, 100, 100); margin: 20px 0 20px 20px\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Invariant zeros</caption>", fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> invariant zero </td>\n  <td> Time constant [s] </td>"
             + "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n</tr>",
            fileName);
        else
          print("<p>\nThe system has no invariant zeros.\n</p>", fileName);
        end if;

      end printHead4;

      encapsulated function printHTMLbasics
        "Print the html preamble or ending on file"
        import Modelica.Utilities.Files;
        import Modelica.Utilities.Streams;

        input String fileName="systemReport.html"
          "File on which the html basics should be written";
        input Boolean printBegin=false
          "True, if beginning of a html file should be printed, otherwise the ending"
          annotation (choices(checkBox=true));

      algorithm
        if printBegin then
          // First print of HTML environment into new file
          Files.removeFile(fileName);
          // Following doesn't work in Dymola
          //Streams.print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", fileName);
          Streams.print("<html>", fileName);
          Streams.print(
            "<head>\n  <title>Analysis of a state space system from Modelica LinearSystems2</title>\n</head>",
            fileName);
          Streams.print("<style type=\"text/css\">", fileName);
          Streams.print("* { font-size: 10pt; font-family: Arial,sans-serif; }",
            fileName);
          Streams.print("</style>", fileName);
        else
          // Last print of HTML environment
          Streams.print("</html>", fileName);
        end if;
      end printHTMLbasics;

      encapsulated function printTab1
        "Print the table with real eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;

        input Eigenvalue evSorted[:];
        input Integer evIndex[size(evSorted, 1)];
        input Real r_evec[size(evSorted, 1), size(evSorted, 1)];
        input Real l_evec[size(evSorted, 1), size(evSorted, 1)];
        input Integer nReal;
        input String xNames2[size(evSorted, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evSorted, 1);
        Real w;
        Real d;

        Integer i;
        Integer j;
        Integer k;
        String number;

        Real r_abs_evec[nx];
        Real l_abs_evec[nx];
        Integer r_maxIndex1;
        Integer l_maxIndex1;
        Integer r_maxIndex2;
        Integer l_maxIndex2;
        //  Complex v_normalized[size(evSorted,1)];
        Real r_abs_v_normalized;
        Real l_abs_v_normalized;
        Real r_v;
        Real l_v;
        Real r_absMax1;
        Real l_absMax1;
        Real r_absMax2;
        Real l_absMax2;
        Boolean r_two;
        Boolean l_two;
        Boolean r_first;
        Boolean l_first;

      algorithm
        i := 1;
        j := i;
        while i <= nReal loop
          // Build eigenvalue number

          number := String(
                  i,
                  minimumLength=7,
                  leftJustified=false);
          j := j + 1;

          // Determine largest value in eigenvector
          k := evIndex[i] "Index with respect to unsorted eigen values";
          r_abs_evec := abs(r_evec[:, k]);
          l_abs_evec := abs(l_evec[:, k]);

          r_first := true;
          r_two := false;
          r_absMax1 := 0;
          r_maxIndex1 := 0;
          r_absMax2 := 0;
          r_maxIndex2 := 0;
          r_abs_v_normalized := Modelica.Math.Vectors.norm(r_abs_evec, 1);

          l_first := true;
          l_two := false;
          l_absMax1 := 0;
          l_maxIndex1 := 0;
          l_absMax2 := 0;
          l_maxIndex2 := 0;
          l_abs_v_normalized := Modelica.Math.Vectors.norm(l_abs_evec, 1);
          for j in 1:nx loop
            r_v := r_abs_evec[j];
            l_v := l_abs_evec[j];

            if r_first then
              r_first := false;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif not r_two then
              r_two := true;
              if r_v < r_absMax1 then
                r_absMax2 := r_v;
                r_maxIndex2 := j;
              else
                r_absMax2 := r_absMax1;
                r_maxIndex2 := r_maxIndex1;
                r_absMax1 := r_v;
                r_maxIndex1 := j;
              end if;
            elseif r_v > r_absMax1 then
              r_absMax2 := r_absMax1;
              r_maxIndex2 := r_maxIndex1;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif r_v > r_absMax2 then
              r_absMax2 := r_v;
              r_maxIndex2 := j;
            end if;

            if l_first then
              l_first := false;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif not l_two then
              l_two := true;
              if l_v < l_absMax1 then
                l_absMax2 := l_v;
                l_maxIndex2 := j;
              else
                l_absMax2 := l_absMax1;
                l_maxIndex2 := l_maxIndex1;
                l_absMax1 := l_v;
                l_maxIndex1 := j;
              end if;
            elseif l_v > l_absMax1 then
              l_absMax2 := l_absMax1;
              l_maxIndex2 := l_maxIndex1;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif l_v > l_absMax2 then
              l_absMax2 := l_v;
              l_maxIndex2 := j;
            end if;

          end for;

          r_absMax1 := 100*r_absMax1/r_abs_v_normalized;
          r_absMax2 := 100*r_absMax2/r_abs_v_normalized;

          if r_absMax2 < 0.05*r_absMax1 then
            r_two := false;
          end if;

          l_absMax1 := 100*l_absMax1/l_abs_v_normalized;
          l_absMax2 := 100*l_absMax2/l_abs_v_normalized;

          if l_absMax2 < 0.05*l_absMax1 then
            l_two := false;
          end if;

          // Print data for one eigen value
          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:center\"> "
             + number + " </td>\n  <td style=\"text-align:left\"> &nbsp; " +
            String(evSorted[i].ev.re, format="14.4e") +
            " </td>\n  <td style=\"text-align:left\"> &nbsp; " + (if evSorted[i].timeConstant
             < 1e6 then String(evSorted[i].timeConstant, format="9.4f") else
            "---") + " </td>\n  <td style=\"text-align:left\"> &nbsp; " + (if
            evSorted[i].isStable then "" else "not ") + "stable, " + (if
            evSorted[i].isStable then (if evSorted[i].isControllable then ""
             else "not ") + "controllable, " else (if evSorted[i].isStabilizable
             then "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable
             then (if evSorted[i].isObservable then "" else "not ") +
            "observable " else (if evSorted[i].isDetectable then "" else "not ")
             + "detectable ") + " </td>", fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td style=\"text-align:left\"> &nbsp; " + " z[" + String(i)
               + "]" + " contributes to " + xNames2[r_maxIndex1] + " with " +
              String(r_absMax1, format=".3g") + " %<br>" + (if r_two then
              "&nbsp; " + " z[" + String(i) + "]" + " contributes to " +
              xNames2[r_maxIndex2] + " with " + String(r_absMax2, format=".3g")
               + " %" else "") + " </td>", fileName);
          end if;

          print("</tr>", fileName);

          i := j;
        end while;

        print("</table>", fileName);
      end printTab1;

      encapsulated function printTab2
        "Print the table with complex conjugate eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;

        input Eigenvalue evSorted[:];
        input Integer evIndex[size(evSorted, 1)];
        input Real r_evec[size(evSorted, 1), size(evSorted, 1)];
        input Real l_evec[size(evSorted, 1), size(evSorted, 1)];
        input Integer nReal;
        input String xNames2[size(evSorted, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evSorted, 1);

        Integer i;
        Integer k;
        String number;
        String number2;
        Integer j;

        Real r_abs_evec[nx];
        Real l_abs_evec[nx];
        Integer r_maxIndex1;
        Integer r_maxIndex2;
        Integer l_maxIndex1;
        Integer l_maxIndex2;
        Real r_abs_v_normalized;
        Real l_abs_v_normalized;
        Real r_v;
        Real l_v;
        Real r_absMax1;
        Real r_absMax2;
        Real l_absMax1;
        Real l_absMax2;
        Boolean r_two;
        Boolean l_two;
        Boolean r_first;
        Boolean l_first;

      algorithm
        i := nReal + 1;
        j := i;
        while i <= nx loop
          // Build eigenvalue number
          number := String(i) + "/" + String(i + 1);
          number2 := number;
          number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;
          j := j + 2;

          // Determine largest value in eigenvector
          k := evIndex[i] "Index with respect to unsorted eigen values";

          for i2 in 1:nx loop
            r_abs_evec[i2] := sqrt(r_evec[i2, k]^2 + r_evec[i2, k + 1]^2);
            l_abs_evec[i2] := sqrt(l_evec[i2, k]^2 + l_evec[i2, k + 1]^2);
          end for;

          r_first := true;
          r_two := false;
          r_absMax1 := 0;
          r_maxIndex1 := 0;
          r_absMax2 := 0;
          r_maxIndex2 := 0;
          r_abs_v_normalized := Modelica.Math.Vectors.norm(r_abs_evec, 1);
          l_first := true;
          l_two := false;
          l_absMax1 := 0;
          l_maxIndex1 := 0;
          l_absMax2 := 0;
          l_maxIndex2 := 0;
          l_abs_v_normalized := Modelica.Math.Vectors.norm(l_abs_evec, 1);

          for j in 1:nx loop
            r_v := r_abs_evec[j];
            l_v := l_abs_evec[j];

            if r_first then
              r_first := false;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif not r_two then
              r_two := true;
              if r_v < r_absMax1 then
                r_absMax2 := r_v;
                r_maxIndex2 := j;
              else
                r_absMax2 := r_absMax1;
                r_maxIndex2 := r_maxIndex1;
                r_absMax1 := r_v;
                r_maxIndex1 := j;
              end if;
            elseif r_v > r_absMax1 then
              r_absMax2 := r_absMax1;
              r_maxIndex2 := r_maxIndex1;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif r_v > r_absMax2 then
              r_absMax2 := r_v;
              r_maxIndex2 := j;
            end if;

            if l_first then
              l_first := false;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif not l_two then
              l_two := true;
              if l_v < l_absMax1 then
                l_absMax2 := l_v;
                l_maxIndex2 := j;
              else
                l_absMax2 := l_absMax1;
                l_maxIndex2 := l_maxIndex1;
                l_absMax1 := l_v;
                l_maxIndex1 := j;
              end if;
            elseif l_v > l_absMax1 then
              l_absMax2 := l_absMax1;
              l_maxIndex2 := l_maxIndex1;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif l_v > l_absMax2 then
              l_absMax2 := l_v;
              l_maxIndex2 := j;
            end if;

          end for;
          r_absMax1 := 100*r_absMax1/r_abs_v_normalized;
          r_absMax2 := 100*r_absMax2/r_abs_v_normalized;
          if r_absMax2 < 0.05*r_absMax1 then
            r_two := false;
          end if;

          l_absMax1 := 100*l_absMax1/l_abs_v_normalized;
          l_absMax2 := 100*l_absMax2/l_abs_v_normalized;
          if l_absMax2 < 0.05*l_absMax1 then
            l_two := false;
          end if;

          // Print data for one eigen value
          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> "
             + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " &plusmn; " +
            String(evSorted[i].ev.im, format="12.4e") + "j" + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(evSorted[i].frequency,
            format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(evSorted[i].damping,
            format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + (if evSorted[i].isStable
             then "" else "not ") + "stable, " + (if evSorted[i].isStable then
            (if evSorted[i].isControllable then "" else "not ") +
            "controllable, " else (if evSorted[i].isStabilizable then "" else
            "not ") + "stabilizable, ") + (if evSorted[i].isStable then (if
            evSorted[i].isObservable then "" else "not ") + "observable " else
            (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td>", fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td style=\"text-align:left\"> &nbsp; " + " z[" + number2
               + "]" + " contribute to " + xNames2[r_maxIndex1] + " with " +
              String(r_absMax1, format=".3g") + " %<br>" + (if r_two then
              "&nbsp; " + " z[" + number2 + "]" + " contribute to " + xNames2[
              r_maxIndex2] + " with " + String(r_absMax2, format=".3g") + " %"
               else "") + " </td>", fileName);
          end if;

          print("</tr>", fileName);

          i := j;
        end while;

        print("</table>", fileName);
      end printTab2;

      encapsulated function printTab3
        "Print the table with eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica.Utilities.Streams.print;
        import Complex;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import LS2Complex = Modelica_LinearSystems2.Math.ComplexAdvanced;

        input Eigenvalue evSorted[:];
        input Complex evecComplex[:, :];
        input Integer evIndex[size(evecComplex, 1)];
        input Complex cev[size(evecComplex, 1)];
        input Integer nReal;
        input String xNames2[size(evecComplex, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evecComplex, 1);

        Integer maxIndex1;
        Integer maxIndex2;

        Complex v_normalized[size(evecComplex, 1)];
        Real abs_v_normalized;
        Real v;
        Real absMax1;
        Real absMax2;
        Boolean two;
        Boolean first;
        Integer j;
        Integer k;
        Integer iw1;
        Integer iw2;
        String number1;
        String number2;
        Real w1;
        Real w2;
        Real d1;
        Real d2;

      algorithm
        for i in 1:nx loop
          // Normalize i-th row of complex eigenvector matrix and determine two largest elements
          v_normalized :=Modelica.ComplexMath.Vectors.normalize(evecComplex[i, :]);
          first := true;
          two := false;
          absMax1 := 0;
          maxIndex1 := 0;
          absMax2 := 0;
          maxIndex2 := 0;
          j := 1;
          abs_v_normalized :=Modelica.ComplexMath.Vectors.norm(v_normalized, 1);
          while j <= nx loop
            if cev[j].im == 0 then
              v := abs(v_normalized[j].re);
              k := j;
              j := j + 1;
            else
              v :=2*Modelica.ComplexMath.'abs'(v_normalized[j]);
              k := j;
              j := j + 2;
            end if;

            if first then
              first := false;
              absMax1 := v;
              maxIndex1 := k;
            elseif not two then
              two := true;
              if v < absMax1 then
                absMax2 := v;
                maxIndex2 := k;
              else
                absMax2 := absMax1;
                maxIndex2 := maxIndex1;
                absMax1 := v;
                maxIndex1 := k;
              end if;
            elseif v > absMax1 then
              absMax2 := absMax1;
              maxIndex2 := maxIndex1;
              absMax1 := v;
              maxIndex1 := k;
            elseif v > absMax2 then
              absMax2 := v;
              maxIndex2 := k;
            end if;
          end while;

          if abs_v_normalized > 1e-30 then
            absMax1 := absMax1/abs_v_normalized;
            absMax2 := absMax2/abs_v_normalized;
          end if;

          if absMax2 < 0.05*absMax1 then
            two := false;
          end if;

          // Determine frequency and number of corresponding eigenvalue
          (w1,d1) := LS2Complex.frequency(cev[maxIndex1]);
          iw1 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex1, evIndex);
          if iw1 <= nReal then
            number1 := String(iw1);
          else
            number1 := String(iw1) + "/" + String(iw1 + 1);
          end if;

          if two then
            (w2,d2) := LS2Complex.frequency(cev[maxIndex2]);
            iw2 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex2, evIndex);
            if iw2 <= nReal then
              number2 := String(iw2);
            else
              number2 := String(iw2) + "/" + String(iw2 + 1);
            end if;
          end if;

          if two then
            print(
              "<tr style=\"background-color:white\">\n  <td rowspan=2 style=\"text-align:left\"> &nbsp; "
               + xNames2[i] + " </td>" +
              "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax1, format="5.1f") + "% by z[" + number1 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number1 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(w1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(d1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then String(evSorted[i].timeConstant, format="9.4f") else
              "--- </td>") + "\n</tr>\n<tr style=\"background-color:white\">"
               + "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax2, format="5.1f") + "% by z[" + number2 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number2 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw2 <= nReal
               then "---" else String(w2, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw2 <= nReal
               then "---" else String(d2, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if (iw2 <=
              nReal and abs(cev[maxIndex2].re) > 1e-10) then String(1/abs(cev[
              maxIndex2].re), format="9.4f") else "--- </td>\n</tr>"), fileName);
          else
            print(
              "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
               + xNames2[i] + " </td>" +
              "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax1, format="5.1f") + "% by z[" + number1 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number1 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(w1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(d1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then String(evSorted[i].timeConstant, format="9.4f") else
              "--- </td>\n</tr>"), fileName);
          end if;
          //     print("<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; " + xNames2[i] + " </td>\n  <td style=\"text-align:left\"> &nbsp; "
          //        + " is composed of " + String(100*absMax1, format="5.1f") + "% by z[" +
          //       number1 + "]" + (if two then " <br>" + " &nbsp; " + " is composed of " +
          //       String(100*absMax2, format="5.1f") + "% by z[" + number2 + "]" else "") + " </td> <td style=\"text-align:center\"> &nbsp; "
          //        + number1 + (if two then "<br> &nbsp; " + number2 else Strings.repeat(9))
          //        + " </td> <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal then
          //             "---" else String(w1, format="9.4f")) + (if two then "<br> &nbsp; "
          //        + (if iw2 <= nReal then "---" else String(w2, format="9.4f")) else
          //       Strings.repeat(9)) + " </td>\n  <td style=\"text-align:center\"> &nbsp; " +
          //       (if iw1 <= nReal then "---" else String(d1, format="9.4f")) + (if two then
          //             "<br> &nbsp; " + (if iw2 <= nReal then "---" else String(d2,
          //       format="9.4f")) else "") + " </td>\n  <td style=\"text-align:center\"> &nbsp; "
          //        + (if (iw1 <= nReal) then String(evSorted[i].timeConstant, format="9.4f") else
          //             "---") + (if two then "<br> &nbsp; " + (if (iw2 <= nReal and abs(
          //       cev[maxIndex2].re) > 1e-10) then String(1/abs(cev[maxIndex2].re),
          //       format="9.4f") else "---") else "") + " </td>\n</tr> ", fileName);

        end for;
        print("</table>", fileName);

      end printTab3;

      encapsulated function printTab4
        "Print the table with eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica.Utilities.Streams.print;
        import Complex;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import LS2Complex = Modelica_LinearSystems2.Math.ComplexAdvanced;

        input Complex systemZeros[:];
        input Integer evIndex[size(systemZeros, 1)];
        input Integer nReal;
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nz=size(systemZeros, 1);

        String number;
        Real timeConstant;
        Real freq;
        Real damp;

      algorithm
        for i in 1:nReal loop
          // Build eigenvalue number

          number := String(
                  i,
                  minimumLength=7,
                  leftJustified=false);
          timeConstant := if abs(systemZeros[i].re) > 10*Modelica.Constants.eps
             then 1/abs(systemZeros[i].re) else 1/(10*Modelica.Constants.eps);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td> &nbsp; " + String(systemZeros[i].re,
            format="14.4e") + " </td>" + "\n  <td> &nbsp; " + String(
            timeConstant, format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>\n</tr>",
            fileName);

        end for;

        for i in nReal + 1:2:nz loop
          number := String(i) + "/" + String(i + 1);
          number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;

          // Determine frequency and number of corresponding zero
          (freq,damp) := LS2Complex.frequency(systemZeros[i]);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(systemZeros[i].re, format="14.4e") + " &plusmn; " +
            String(systemZeros[i].im, format="12.4e") + "j </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(freq, format=
            "9.4f") + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(damp, format="9.4f") + " </td>\n</tr>", fileName);

        end for;

        print("</table>\n", fileName);
      end printTab4;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(ss);
   or
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(
  ss,
  analyseOptions=<a href=\"modelica://Modelica_LinearSystems2.Internal.AnalyseOptions\">analyseOptions</a>,
  fileName,
  systemName,
  description);
</pre></blockquote>

<h4>Description</h4>
<p>
This function analyzes a state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>     <label for=\"eqn1\">(1)</label>
    <b>x</b>(t=0) = <b>x</b><sub>0</sub>
</pre></blockquote>
<p>
based on its poles, i.e. the eigenvalues, and the zeros of the system.
The system will be checked for stability, controllability and observability. In the case that the system is not stable stabilizability and detectability are examined. Furthermore, stability, controllability, observability, stabilizability, and detectability are indicated for each eigenvalue.
</p>

<h5>Stability</h5>
<p>
System (1) is stable if and only if all eigenvalues of the matrix <b>A</b> have negative real parts.
The calculation of the eigenvalues is based on the LAPACK routine dgeev.
</p>

<h5>Controllability</h5>
<p>
System (1) is said to be controllable if, starting from any initial state <b>x</b><sub>0</sub>, the system can be driven by appropriate inputs to any final state <b>x</b><sub>1</sub> within some finite time window. Equivalent is that the eigenvalues of <b>A</b>-<b>BK</b> can  arbitrarily be assigned by an appropriate choice of the matrix <b>K</b>.
</p>

<h5>Stabilizability</h5>
<p>
System (1) is said to be stabilizable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are controllable. Therefore, a controllable system is always stabilizable. An equivalent definition of stabilizability is, that a system is said to be stabilizable if there exist a matrix <b>K</b> such that <b>A</b>-<b>BK</b> is stable.
</p>

<h5>Observability</h5>
<p>
System (1) is said to be observable if the (arbitrary) initial state <b>x</b><sub>0</sub> can be uniquely determined from any state <b>x</b>(t<sub>1</sub>), t<sub>1</sub>>0, from the knowledge of the input <b>u</b>(t) and output <b>y</b>(t). With other words,  from the system's outputs it is possible to determine the behavior of the entire system. Equivalent is, that the eigenvalues of <b>A</b>-<b>LC</b> can be arbitrarily be assigned by an appropriate choice of matrix <b>L</b>.
Observability is called the dual concept of controllability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is observable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is controllable.
</p>

<h5>Detectability</h5>
<p>
System (1) is said to be detectable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are observable. Therefore, a observable system is always detectable. An equivalent definition of detectability is, that a system is said to be detectable if there exist a matrix <b>L</b> such that <b>A</b>-<b>LC</b> is stable.
Detectability is called the dual concept of stabilizability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is detectable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is stabilizable.
</p>

<h5>Algorithm to test controllability/stabilizability and observability/detectability respectively</h5>
<p>
The test of controllability and stabilizability is performed with the staircase algorithm which transforms the system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) into the controller-Hessenberg form (<b>A</b><sub>H</sub>, <b>B</b><sub>H</sub>, <b>C</b><sub>H</sub>, <b>D</b>) with <b>A</b><sub>H</sub> is a block upper Hessenberg matrix and <b>B</b><sub>H</sub>=[<b>B</b><sub>1</sub>; 0] with triangular matrix <b>B</b><sub>1</sub> with rank(<b>B</b><sub>1</sub>) = rank(<b>B</b>).
In <b>A</b><sub>H</sub>=[<b>A</b><sub>c</sub>, *,0, <b>A</b><sub>nc</sub>) the eigenvalues of the matrices <b>A</b><sub>c</sub> and <b>A</b><sub>nc</sub> are the controllable eigenvalues and uncontrollable eigenvalues of <b>A</b> respectively.
The test of observability and detectability is performed by testing the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) with respect to controllability and stabilizability.
</p>

<h5>Solution of a linear time invariant system </h5>
<p>
The solution <b>x</b>(t) of the initial value problem (1) consists of the homogeneous part (zero input response) <b>x</b><sub>h</sub>(t) and the inhomogeneous part x<sub>i</sub>(t). The zero input solution is given by
</p>
<blockquote><pre>
<b>x</b><sub>h</sub>(t) = exp(<b>A</b>*(t-t<sub>0</sub>))<b>x</b><sub>0</sub>.
</pre></blockquote>
<p>
The system can also be represented as a linear combination of the modal states <b>z</b>,
</p>
<blockquote><pre>
<b>x</b> = <b>V</b><b>z</b>
</pre></blockquote>
<p>
i.e. the states of a similar system, with
</p>
<blockquote><pre>
der(<b>z</b>) = <b>V</b><sup>-1</sup><b>AVz</b> + <b>V</b><sup>-1</sup><b>B</b><b>u</b>
</pre></blockquote>
<p>
where the system matrix <b>V</b><sup>-1</sup><b>AV</b> is the real Jordan form. For single real eigenvectors the system is decoupled, i.e. the solution of the modal states are denoted by
<blockquote><pre>
z<sub>i</sub> = exp(s<sub>i</sub> t)*z<sub>0i</sub>
</pre></blockquote>
<p>
The behavior of the modal states is determined as the solution of a linear first order differential equation for real eigenvalues. Since this behavior is well known, the behavior of the x<sub>i</sub> can at least roughly be estimated by means of the behavior of the most relevant modal states. Therefore, the contribution of the modal states to the states is computed as an indication of the original system behavior.
</p>

<h5>Contribution of the modal states to the states</h5>
<p>
Generally, as described above, the states of the system can be described as linear combination of modal states and, therefore, the states can be characterized to a certain extend by the modal states if the proportions of the combination are known. Hence, for each modal state z<sub>i</sub> of the vector <b>z</b> the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub>| of the corresponding right eigenvector <b>v</b><sub>i</sub> indicate the proportion of <b>z</b><sub>i</sub> that is contributed to the state x<sub>j</sub>.
On the other hand, the composition of xi is indicated by the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>|, i.e. the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>| of the corresponding row <b>v</b><sub>i</sub><sup>T</sup> of the eigenvector matrix <b>V</b> indicate the proportion of the state x<sub>i</sub> that is contributed by the modal state z<sub>j</sub>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ss=StateSpace(
    A=[-3,2,-3,4,5,6; 0,6,7,8,9,4; 0,2,3,0,78,6; 0,1,2,2,3,3; 0,13,34,0,0,1; 0,
      0,0,-17,0,0],
    B=[1,0; 0,1; 1,0; 0,1; 1,0; 0,1],
    C=[0,0,1,0,1,0; 0,1,0,0,1,1],
    D=[0,0; 0,0],
    xNames={\"x1\",\"x2\",\"x3\",\"x4\",\"x5\",\"x6\"},
    uNames={\"u1\",\"u2\"}, yNames={\"y1\",\"y2\"});

  String fileName=\"analysis.html\";
  String systemName=\"Demonstration System\";
  String description=\"System to demonstrate the usage of Modelica_LinearSystems2.StateSpace.Analysis.anlysis()\"

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Analysis.analysis(ss, fileName=fileName, systemName=systemName, description=description)
//  gives:
</pre></blockquote>

<h4>System report</h4>
<p>
The system <b>Demonstation System</b>
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>
</pre></blockquote>
<p>
is defined by
</p>
<blockquote><pre>
        x1   x2   x3   x4   x5   x6            u1  u2
    x1  -3    2   -3    4    5    6         x1  1   0
    x2   0    6    7    8    9    4         x2  0   1
A = x3   0    2    3    0   78    6     B = x3  1   0
    x4   0    1    2    2    3    3         x4  0   1
    x5   0   13   34    0    0    1         x5  1   0
    x6   0    0    0  -17    0    0         x6  0   1

        x1   x2   x3   x4   x5   x6            u1  u2
C = y1   0    0    1    0    1    0     D = y1  0   0
    y2   0    1    0    0    1    1         y2  0   0
</pre></blockquote>

<h5>Description</h5>
<p>
System to demonstrate the usage of Modelica_LinearSystems2.StateSpace.Analysis.analysis()
</p>

<h5>Characteristics</h5>
<p>The system
<br> is
not stable
<br>but it is controllable
<br> and therefore it is stabilizable
<br> The system is not observable
<br> but it is detectable
</p>

<p>
<b><big>Eigenvalues analysis</big></b>
<br><br>
<b>Real eigenvalues</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin: 20px 0 20px 20px\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td><td> eigenvalue </td> <td> T [s] </td>  <td> characteristics </td><td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:center\">       1 </td> <td style=\"text-align:left\"> &nbsp;   -4.9874e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0201 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, observable  </td> <td style=\"text-align:left\"> &nbsp;  z[1] contributes to x3 with 54.6 %<br>&nbsp;  z[1] contributes to x5 with 37 % </td> </tr>
<tr>
 <td style=\"text-align:center\">       2 </td> <td style=\"text-align:left\"> &nbsp;   -3.0000e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3333 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, not observable  </td> <td style=\"text-align:left\"> &nbsp;  z[2] contributes to x1 with 100 %<br> </td> </tr>
<tr>
 <td style=\"text-align:center\">       3 </td> <td style=\"text-align:left\"> &nbsp;    2.9891e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3346 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[3] contributes to x2 with 51.9 %<br>&nbsp;  z[3] contributes to x1 with 23.9 % </td> </tr>
<tr>
 <td style=\"text-align:center\">       4 </td> <td style=\"text-align:left\"> &nbsp;    5.5825e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0179 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[4] contributes to x3 with 48.4 %<br>&nbsp;  z[4] contributes to x5 with 32.5 % </td> </tr>
</table>

<p>
<b>Conjugated complex pairs of eigenvalues</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin: 20px 0 20px 20px\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> eigenvalue </td><td> freq. [Hz] </td> <td> damping </td><td> characteristics </td>  <td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:left\">     5/6 </td> <td style=\"text-align:left\"> &nbsp;    1.0299e+000 &plusmn;  6.5528e+000j </td> <td style=\"text-align:left\"> &nbsp;    1.0557 </td> <td style=\"text-align:left\"> &nbsp;   -0.1553 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[    5/6] contribute to x6 with 35.9 %<br>&nbsp;  z[    5/6] contribute to x2 with 20.6 % </td> </tr>
</table>

<p>
In the table above, the column <b>contribution to states</b> lists for each eigenvalue the states
to which thecorresponding modal state contributes most. This information is based on the
two largest absolute values of the corresponding right eigenvector (if the second large value
is less than 5&nbsp;% of the largest contribution, it is not shown).
</p>

<p>
In the next table, for each state in the column <b>correlation to modal states</b>, the modal
states which contribute most to the coresponding state are summarized, i.e. the state is mostly composed of these modal states
This information is based on the two largest absolute values of row i of the
eigenvector matrix that is associated with eigenvalue i (if the second large value
is less than 5&nbsp;% of the largest contribution, it is not shown). This only holds
if the modal states are in the same order of magnitude. Otherwise, the modal states
listed in the last column might be not the most relevant one.
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin: 20px 0 20px 20px\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> state </td> <td> composition </td> <td> eigenvalue #</td> <td> freq. [Hz] </td> <td> damping </td> <td> T [s] </td></tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x1 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  42.5% by z[2] <br> &nbsp;  is composed of  35.4% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 2<br> &nbsp; 5/6 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    1.0557 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp;    0.0201<br> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x2 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  44.2% by z[3] <br> &nbsp;  is composed of  43.7% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 3<br> &nbsp; 5/6 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    1.0557 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp;    0.3333<br> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x3 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  36.9% by z[1] <br> &nbsp;  is composed of  36.3% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 1<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;    0.3346<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x4 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  88.9% by z[5/6] <br> &nbsp;  is composed of   9.8% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 5/6<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp;    1.0557<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;   -0.1553<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x5 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  45.3% by z[1] <br> &nbsp;  is composed of  44.1% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 1<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;    0.0000<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x6 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  95.7% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 5/6          </td> <td style=\"text-align:center\"> &nbsp;    1.0557          </td> <td style=\"text-align:center\"> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
</table>

<p>
<b>Invariant zeros</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin: 20px 0 20px 20px\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> invariant zero </td><td> Time constant [s] </td> <td> freq. [Hz] </td> <td> damping </td></tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;       1 </td> <td> &nbsp;   -5.4983e+001 </td> <td> &nbsp;    0.0182 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;       2 </td> <td> &nbsp;   -3.0000e+000 </td> <td> &nbsp;    0.3333 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;     3/4 </td> <td style=\"text-align:left\"> &nbsp;    3.2417e+000 &plusmn;  5.6548e+000j </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:left\"> &nbsp;    1.0374 </td> <td style=\"text-align:left\"> &nbsp;   -0.4973 </td> </tr>
</table>
</html>", revisions="<html>
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
      "Calculate the time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      input TimeResponse response=TimeResponse.Step;
      extends Modelica_LinearSystems2.Internal.timeResponseMask2(redeclare Real
          y[:, size(sc.C, 1), if response == TimeResponse.Initial then 1 else
          size(sc.B, 2)], redeclare Real x_continuous[:, size(sc.A, 1), if
          response == TimeResponse.Initial then 1 else size(sc.B, 2)]);
      // Input/Output declarations of time response functions

      input Real x0[size(sc.A, 1)]=zeros(size(sc.A, 1)) "Initial state vector";

    protected
      Real dtVar;
      Real tSpanVar;
      Integer samples;
      Real u[:, size(sc.B, 2)];
      Real new_x[size(sc.A, 1), 1];
      Real x[size(sc.A, 1), 1]=zeros(size(sc.A, 1), 1);
      Modelica_LinearSystems2.DiscreteStateSpace sd(
        redeclare Real A[size(sc.A, 1), size(sc.A, 2)],
        redeclare Real B[size(sc.B, 1), size(sc.B, 2)],
        redeclare Real C[size(sc.C, 1), size(sc.C, 2)],
        redeclare Real D[size(sc.D, 1), size(sc.D, 2)],
        redeclare Real B2[size(sc.B, 1), size(sc.B, 2)]);
      Real i1;
      Real i2;

    algorithm
      // Check whether system has inputs and outputs
      assert(size(sc.B, 2) > 0,
        "\n... StateSpace system has no inputs and it is not possible to compute a time response.\n");
      assert(size(sc.C, 1) > 0,
        "\n... StateSpace system has no outputs and it is not possible to compute a time response.\n");

      // set sample time and simulation time span
      if (dt == 0 and tSpan == 0) then
        (dtVar,tSpanVar) :=
          Modelica_LinearSystems2.Internal.timeResponseSamples(sc);
      elseif (dt == 0 and tSpan <> 0) then
        dtVar := Modelica_LinearSystems2.Internal.timeResponseSamples(sc);
        tSpanVar := tSpan;
      elseif (dt <> 0 and tSpan == 0) then
        (,tSpanVar) := Modelica_LinearSystems2.Internal.timeResponseSamples(sc);
        dtVar := dt;
      else
        dtVar := dt;
        tSpanVar := tSpan;
      end if;

      samples := integer(tSpanVar/dtVar + 1);
      t := 0:dtVar:tSpanVar;
      u := zeros(samples, size(sc.B, 2));
      y := if response == TimeResponse.Initial then zeros(
            samples,
            size(sc.C, 1),
            1) else zeros(
            samples,
            size(sc.C, 1),
            size(sc.B, 2));
      x_continuous := if response == TimeResponse.Initial then zeros(
            samples,
            size(sc.A, 1),
            1) else zeros(
            samples,
            size(sc.A, 1),
            size(sc.B, 2));

      if response == TimeResponse.Initial then
        sd :=Modelica_LinearSystems2.DiscreteStateSpace(
              sc,
              dtVar,
              Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal);
        (y[:, :, 1],x_continuous[:, :, 1]) :=
          Modelica_LinearSystems2.DiscreteStateSpace.initialResponse(
              sd,
              x0,
              samples);
      else

        for i1 in 1:size(sc.B, 2) loop
          // Loop over inputs

          // time response to plot
          if response == TimeResponse.Impulse then
            u[1, :] := zeros(size(sc.B, 2));
            u[1, i1] := 1;
            sd :=Modelica_LinearSystems2.DiscreteStateSpace(
                  sc,
                  dtVar,
                  Modelica_LinearSystems2.Utilities.Types.Method.ImpulseExact);
          elseif response == TimeResponse.Step then
            u[:, :] := zeros(samples, size(sc.B, 2));
            u[:, i1] := ones(samples);
            sd :=Modelica_LinearSystems2.DiscreteStateSpace(
                  sc,
                  dtVar,
                  Modelica_LinearSystems2.Utilities.Types.Method.StepExact);
          elseif response == TimeResponse.Ramp then
            u[:, :] := zeros(samples, size(sc.B, 2));
            u[:, i1] := 0:dtVar:tSpanVar;
            sd :=Modelica_LinearSystems2.DiscreteStateSpace(
                  sc,
                  dtVar,
                  Modelica_LinearSystems2.Utilities.Types.Method.RampExact);
            //    elseif response == TimeResponse.Initial then
            //      u[:, :] := zeros(samples, size(sc.B, 2));
            //      sd := Modelica_LinearSystems2.DiscreteStateSpace(
            //          sc,
            //          dtVar,
            //          Modelica_LinearSystems2.Types.Method.Trapezoidal);
          else
            assert(false, "Argument response (= " + String(response) +
              ") of \"Time response to plot\" is wrong.");
          end if;

          (y[:, :, i1],x_continuous[:, :, i1]) :=
            Modelica_LinearSystems2.DiscreteStateSpace.timeResponse(
                sd,
                u,
                x0);
        end for;
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = StateSpace.Analysis.<b>timeResponse</b>(ss, dt, tSpan, responseType, x0)
</pre> </blockquote>

<h4>Description</h4>
<p>
This function calculates the time responses of a state space system. The type of the time response is defined by the input <code>responseType</code>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step &quot;Step response&quot;,
Ramp &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;
</pre></blockquote>
<p>
The state space system is transformed to a appropriate discrete state space system and,
starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and x are calculated for each time step t=k*dt.
</p>

<h4>Example</h4>
<blockquote><pre>  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[2],
    D=[0]);
  Real Ts = 0.1;
  Real tSpan = 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response = Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[1] = {0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(ss,Ts,tSpan,response,x0);
  // y[:,1,1] = {0, 0.19, 0.3625, 0.518, 0.659}
  //        t = {0, 0.1, 0.2, 0.3, 0.4}
  // x[:,1,1] = {0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>
</html>", revisions="<html>
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
    end timeResponse;

    encapsulated function impulseResponse
      "Calculate the impulse time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask2;

    algorithm
      (y,t,x_continuous) :=Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
            sc=sc,
            dt=dt,
            tSpan=tSpan,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
            x0=zeros(size(sc.A, 1)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = StateSpace.Analysis.<b>impulseResponse</b>(ss, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the time response of a state space system for impulse imput. The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
StateSpace.Analysis.impulseResponse(ss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,size(ss.A,1))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[2],
    D=[0]);
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x) := StateSpace.Analysis.impulseResponse(ss,Ts,tSpan);
  // y[:,1,1] = {2, 1.8097, 1.6375, 1.4816, 1.3406}
  //        t = {0, 0.1, 0.2, 0.3, 0.4}
  // x[:,1,1] = {1, 0.9048, 0.8187, 0.7408, 0.6703}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>
</html>", revisions="<html>
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
    end impulseResponse;

    encapsulated function stepResponse
      "Calculate the step time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask2;

    algorithm
      (y,t,x_continuous) :=Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
            sc=sc,
            dt=dt,
            tSpan=tSpan,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
            x0=zeros(size(sc.A, 1)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = StateSpace.Analysis.<b>stepResponse</b>(ss, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>stepResponse</b> calculates the step response of a state space system.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
StateSpace.Analysis.stepResponse(ss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Step, x0=fill(0,size(ss.A,1))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[2],
    D=[0]);
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=StateSpace.Analysis.stepResponse(ss,Ts,tSpan);
//  y[:,1,1]={0, 0.19, 0.3625, 0.518, 0.659}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>
</html>", revisions="<html>
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
    end stepResponse;

    encapsulated function rampResponse
      "Calculate the ramp time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask2;

    algorithm
      (y,t,x_continuous) :=Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
            sc=sc,
            dt=dt,
            tSpan=tSpan,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
            x0=zeros(size(sc.A, 1)));

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = StateSpace.Analysis.<b>rampResponse</b>(ss, dt, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rampResponse</b> calculates the time response of a state space system for ramp imput u = t.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
StateSpace.Analysis.rampResponse(ss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,size(ss.A,1))).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1],
      B=[1],
      C=[2],
      D=[0]);
  Real Ts=0.1;
  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=StateSpace.Analysis.rampResponse(ss,Ts,tSpan);
//  y[:,1,1]={0, 0.00967, 0.03746, 0.08164, 0.14064}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.00484, 0.018734, 0.04082, 0.07032}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>
</html>", revisions="<html>
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
    end rampResponse;

    encapsulated function initialResponse
      "Calculate the time response of a state space system for given initial condition and zero inputs"

      import Modelica;
      import Modelica_LinearSystems2;

      input Real x0[:]=fill(0, 0) "Initial state vector";

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask2(redeclare Real
          y[:, size(sc.C, 1), 1], redeclare Real x_continuous[:, size(sc.A, 1),
          1]);

    algorithm
      (y,t,x_continuous) :=Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
            sc=sc,
            dt=dt,
            tSpan=tSpan,
            response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
            x0=x0);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(y, t, x) = StateSpace.Analysis.<b>initialResponse</b>(ss, dt, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>initialResponse</b> calculates the time response of a state space system for given initial condition and zero inputs.
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
</p>
<blockquote><pre>
StateSpace.Analysis.initialResponse(x0,ss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1],
      B=[1],
      C=[2],
      D=[0]);
  Real Ts=0.1;
  Real tSpan= 0.4;
  Real x0[2] = {1};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=StateSpace.Analysis.initialResponse(x0,ss,Ts,tSpan);
//  y[:,1,1]={2, 1.809, 1.637, 1.4812, 1.3402}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 0.9048, 0.8186, 0.7406, 0.6701}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>
</html>", revisions="<html>
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
    end initialResponse;

    encapsulated function numeratorDegree
      "Return numerator degree of the corresponding transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input Modelica_LinearSystems2.StateSpace ss;
      output Integer result;

    protected
      Modelica_LinearSystems2.TransferFunction tf=if StateSpace.Internal.isSISO(
          ss) then StateSpace.Conversion.toTransferFunction(ss) else
          TransferFunction(1);

    algorithm
      assert(StateSpace.Internal.isSISO(ss), "System must be SISO but is " +
        String(size(ss.B, 2)) + "-by-" + String(size(ss.C, 1)) + " system");
      result := size(tf.n, 1) - 1;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = StateSpace.Analysis.<b>numeratorDegree</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>numeratorDegree</b> calculates the degree of the numerator polynomial of the corresponding transfer function.
The state space system is converted to the transfer function G(s)=N(s)/D(s) with the polynomial N(s) as numerator.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[1],
    D=[1]);

  Real nDegree;

<b>algorithm</b>
  nDegree := StateSpace.Analysis.numeratorDegree(Modelica_LinearSystems2.StateSpace(ss);
//  nDegree = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.denominatorDegree\">StateSpace.Analysis.denominatorDegree</a>
</p>
</html>"));
    end numeratorDegree;

    encapsulated function denominatorDegree
      "Return denominator degree of the corresponding transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;

      input Modelica_LinearSystems2.StateSpace ss;
      output Integer result;
    protected
      Modelica_LinearSystems2.TransferFunction tf=if StateSpace.Internal.isSISO(
          ss) then StateSpace.Conversion.toTransferFunction(ss) else
          TransferFunction(1);

    algorithm
      assert(StateSpace.Internal.isSISO(ss), "System must be SISO but is " +
        String(size(ss.B, 2)) + "-by-" + String(size(ss.C, 1)) + " system");
      result := size(tf.d, 1) - 1;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = StateSpace.Analysis.<b>denominatorDegree</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of the corresponding transfer function.
The state space system is converted to the transfer function G(s)=N(s)/D(s) with the polynomial D(s) as denominator.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[1],
    D=[1]);

  Real dDegree;

<b>algorithm</b>
  dDegree := StateSpace.Analysis.denominatorDegree(Modelica_LinearSystems2.StateSpace(ss);
//  nDegree = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.numeratorDegree\">StateSpace.Analysis.numeratorDegree</a>
</p>
</html>"));
    end denominatorDegree;

    encapsulated function evaluate
      "Evaluate the corresponding transfer function at a given (complex) value of s"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      input Complex s "Value of s where tf shall be evaluated";
      input Real den_min=0 "|denominator(s)| is limited by den_min";
      output Complex result "= tf(s)";

    protected
      Boolean issiso=StateSpace.Internal.isSISO(ss);
      TransferFunction tf=if issiso then
          StateSpace.Conversion.toTransferFunction(ss) else TransferFunction(1);
      Complex j=Modelica.ComplexMath.j;
      Complex den=Polynomial.evaluateComplex(Polynomial(tf.d), s);
      Real abs_den=Modelica.ComplexMath.'abs'(den);
    algorithm
      assert(issiso, "System must be SISO but is " + String(size(ss.B, 2)) +
        "-by-" + String(size(ss.C, 1)) + " system");
      den := if abs_den >= den_min then den else -abs_den + 0*j;
      result := Polynomial.evaluateComplex(Polynomial(tf.n), s)/den;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = StateSpace.Analysis.<b>evaluate</b>(ss,s)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the corresponding transfer function of the state space system at a given (complex) value of s.
The state space system is converted to the transfer function G(s)=N(s)/D(s), which is evaluated by calculating the numerator polynomial N(s) and the denominator polynomial D(s).
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[1],
    D=[0]);
  Complex s=Complex(1,1);

  Complex result;

<b>algorithm</b>
  result := Modelica_LinearSystems2.StateSpace.Analysis.evaluate(ss, s);
//  result = 0.4 - 0.2j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a>,
<a href=\"Modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>
</html>", revisions="<html>
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
      "Calculate zeros and poles of the TransferFunction corresponding to a state space representation"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;

      output Complex z[:] "Zeros (Complex vector of numerator zeros)";
      output Complex p[:] "Poles (Complex vector of denominator zeros)";
      output Real k
        "Constant multiplied with transfer function that is factorized with zeros and poles";

    protected
      TransferFunction tf=StateSpace.Conversion.toTransferFunction(ss);
      Polynomial pn;
      Polynomial pd;
      TransferFunction tf2;
      Real r;
      Complex s;
      Complex y1;
      Complex y2;
    algorithm

      z := Polynomial.roots(Polynomial(tf.n));
      p := Polynomial.roots(Polynomial(tf.d));
      pn := Polynomial(z);
      pd := Polynomial(p);
      tf2 := TransferFunction(pn, pd);
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
(z,p,k) = StateSpace.Analysis.<b>zerosAndPoles</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the zeros, poles and gain of the corresponding transfer function of a state space system.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[1],
    D=[1]);

  Complex z;
  Complex p;
  Real k;

<b>algorithm</b>
  (z,p,k)=Modelica_LinearSystems2.StateSpace.Analysis.zerosAndPoles(ss);
//  z = {-2}
//  p = {-1}
//  k = 1
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a>,
<a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">TransferFunction.Analysis.zerosAndPoles</a>
</p>
</html>"));
    end zerosAndPoles;

    encapsulated function eigenValues
      "Calculate the eigenvalues of a linear state space system and write them in a complex vector"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Complex eigvalues[size(ss.A, 1)]=
        Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A)
        "Eigen values of the system";
    algorithm

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
eigenvalues = StateSpace.Analysis.<b>eigenValues</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvalues of a state space system, i.e. the eigenvalues of the system matrix <b>A</b> of a state space system. The output is a complex vector containing the eigenvalues.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1,1;-1,-1],
    B=[1;1],
    C=[1,1],
    D=[0]);

  Complex eigenvalues[2];

<b>algorithm</b>
  eigenvalues = Modelica_LinearSystems2.StateSpace.Analysis.eigenValues(ss);
// eigenvalues = {-1 + 1j, -1 - 1j}
</pre></blockquote>
</html>"));
    end eigenValues;

    encapsulated function eigenVectors
      "Calculate the rigth eigenvectors of a linear state space system and write them columnwise in a matrix. Optionally, the eigenvalues are computed"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      input Boolean onlyEigenvectors=true "True, if only eigen vertor eigvec should be calculated";
      output Real eigvec[size(ss.A, 1), size(ss.A, 2)]
        "Eigen values of the system";
      output Complex eigval[size(ss.A, 1)]=fill(Complex(0), size(ss.A, 1))
        "Eigen values of the system";
    protected
      Integer info;
      Real eigvalRe[size(ss.A, 1)]=fill(0, size(ss.A, 1));
      Real eigvalIm[size(ss.A, 1)]=fill(0, size(ss.A, 1));

    algorithm
      if size(ss.A, 1) > 0 then
        if onlyEigenvectors then
          (,,eigvec,info) := LAPACK.dgeev(ss.A);
        else
          (eigvalRe,eigvalIm,eigvec,info) := LAPACK.dgeev(ss.A);
          for i in 1:size(ss.A, 1) loop
            eigval[i].re := eigvalRe[i];
            eigval[i].im := eigvalIm[i];
          end for;
        end if;

        assert(info == 0, "Calculating the eigen values with function
\"StateSpace.Analysis.eigenVectors\" is not possible, since the
numerical algorithm does not converge.");
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(eigenvectors, eigenvalues) = StateSpace.Analysis.<b>eigenVectors</b>(ss, onlyEigenvectors)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvectors and optionally (onlyEigenvectors=false) the eigenvalues of a state space system. The output <tt>eigenvectors</tt> is a matrix with the same dimension as matrix <b>ss.A</b>. Just like in <a href=\"modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>, if the i-th eigenvalue has an imaginary part, then <tt>eigenvectors</tt>[:,i] is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the eigenvector of the i-th eigenvalue.<br>
The eigenvalues are returned as a complex vector <tt>eigenvalues</tt>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1,1;-1,-1],
    B=[1;1],
    C=[1,1],
    D=[0]);

  Real eigenvectors[2,2];
  Complex eigenvalues[2];

<b>algorithm</b>
  (eigenvectors, eigenvalues) = Modelica_LinearSystems2.StateSpace.Analysis.eigenVectors(ss, true);
// eigenvectors = [0.707, 0; 0, 0.707]
// eigenvalues = {-1 + 1j, -1 - 1j}

          |0.707 |         | 0.707 |
i.e. v1 = |      |,   v2 = |       |
          |0.707i|         |-0.707i|
</pre></blockquote>
</html>"));
    end eigenVectors;

    encapsulated function invariantZeros
      "Compute invariant zeros of linear state space system"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Complex;

      input StateSpace ss "Linear system in state space form";

      output Complex Zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
      Integer n=10;
      Integer m;
      Integer p;
      Real Ar[:, :];
      Real Br[:, :];
      Real Cr[:, :];
      Real Dr[:, :];

      Real Af[:, :];
      Real Bf[:, :];
      Real AfBf[:, :];

      Real V2[:, :];
      Real Vf[:, :];
      Real R[:, :];

      Integer na;
      Real alphaReal[:];
      Real alphaImag[:];
      Real beta[:];
      Integer info;
      Real zerosMax;
      Real absZero;

      Integer j;
    algorithm
      if min(size(ss.B)) == 0 or min(size(ss.C)) == 0 then
        Zeros := fill(Complex(0), 0);
      else
        (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(
              ss.A,
              ss.B,
              ss.C,
              ss.D);
        if n > 0 then
          (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(
                transpose(Ar),
                transpose(Cr),
                transpose(Br),
                transpose(Dr));
        end if;
        if n == 0 then
          Zeros := fill(Complex(0), 0);
        else
          (,R,,V2) := Matrices.QR(Matrices.fliplr(transpose([Cr, Dr])));
          Vf := Matrices.fliplr(V2);
          AfBf := [Ar, Br]*Vf;
          Af := AfBf[:, 1:size(Ar, 2)];
          Bf := Vf[1:size(Ar, 1), 1:size(Ar, 2)];

          (alphaReal,alphaImag,beta,,,info) := LAPACK.dggev(
                Af,
                Bf,
                n);
          assert(info == 0,
            "Failed to compute invariant zeros with function invariantZeros(..)");

          Zeros := fill(Complex(0), size(beta, 1));

          /* The pencil (Af,Bf) has n eigenvalues, since the transformation to this
             form is done so that Bf is regular. Therefore, the generalized eigenvalues
             represented by alpha[i]/beta[i] have the property that beta[i] cannot be zero
             and a division by beta[i] is uncritical.

             |alpha| / beta <= zerosMax
             if |alpha| <= beta*zerosMax then
                zero is used
             else
                assumed that zero is at infinity (i.e. it is ignored)
          */
          j := 0;
          zerosMax := 1.0e4*Modelica.Math.Matrices.norm([Af, Bf], p=1);
          for i in 1:size(beta, 1) loop
             absZero :=Modelica.ComplexMath.'abs'(Complex(alphaReal[i], alphaImag[i]));
             if absZero <= beta[i]*zerosMax then
                j := j + 1;
                Zeros[j].re := alphaReal[i]/beta[i];
                Zeros[j].im := alphaImag[i]/beta[i];
             end if;
          end for;

          if j == 0 then
             Zeros := fill(Complex(0), 0);
          else
             Zeros := Zeros[1:j];
          end if;
        end if;
      end if;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zeros = StateSpace.Analysis.<b>invariantZeros</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
     <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre></blockquote>
<p>
The invariant zeros of this system are defined as the variables
s  that make the Rosenbrock matrix of the system
</p>
<pre>
    | s<b>I-A</b>   <b>-B</b> |
    |           |
    | <b>C</b>       <b>D</b> |
</pre>
<p>
singular.
</p>
<p>
This function applies the algorithm described in [1] where the system (<b>A</b>, <b>B</b>, <b>C</b>, <b>D</b>) is reduced to a new system (<b>A</b>r, <b>B</b>r <b>C</b>r, <b>D</b>r) with the same zeros and with <b>D</b>r of full rank.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[1, 1, 1;0, 1, 1;0,0,1],
    B=[1;0;1],
    C=[0,1,1],
    D=[0]);

  Complex zeros[:];

<b>algorithm</b>
  zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);
// zeros = {1, 0}
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Emami-Naeini, A. and Van Dooren, P. (1982):</dt>
<dd> <b>Computation of Zeros of Linear Multivariable Systems</b>.
     Automatica, 18, pp. 415-430.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
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
      "Return steady state gain matrix K (for a stable system: K[i,j] = value of y[i] at infinite time for a step input of u[j])"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "Linear system in state space form";

      output Real K[size(ss.C, 1), size(ss.B, 2)]
        "DC gain matrix K (= G(s=0) = D - C*inv(A)*B)";
      output Boolean finite
        "True, if K is finite, otherwise K is infinite (K=fill(Modelica.Constants.inf,..) returned)";

    protected
      Integer nx=size(ss.A, 1);
      Integer nu=size(ss.B, 2);
      Integer ny=size(ss.C, 1);
      Real X[nx, nu];
      Integer rank;
    algorithm
      finite := true;
      if nu == 0 or ny == 0 then
        K := fill(0.0, ny, nu);
      else
        (X,rank) := Modelica_LinearSystems2.Math.Matrices.leastSquares2(ss.A, ss.B);
        // Determine whether A*X-B=0 is not fulfilled (since no unique solution)
        if rank < nx then
          if Modelica.Math.Matrices.norm(ss.A*X - ss.B, p=Modelica.Constants.inf)
               >= 1000*Modelica.Constants.eps then
            finite := false;
          end if;
        end if;

        if finite then
          // A*X - B = 0:
          K := ss.D - ss.C*X;
        else
          // The least squares solution does not fulfill A*X - B = 0
          K := fill(
                Modelica.Constants.inf,
                ny,
                nu);
        end if;
      end if;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
          K = <b>dcGain</b>(ss);
(K, finite) = <b>dcGain</b>(ss);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the steady state gain <b>K</b> of a state space system.
<b>K</b> is defined in the following way:
</p>
<p>
The linear state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
     <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre></blockquote>
<p>
is solved for <b>y</b> under steady state conditions, i.e.,
</p>
<blockquote><pre>
   der(<b>x</b>) = <b>0</b>
</pre></blockquote>
<p>
resulting in
</p>
<blockquote><pre>
    <b>y</b> = ( <b>D</b> + <b>C</b>*inv(<b>A</b>)*<b>B</b> )*<b>u</b>
      = <b>K</b>*<b>u</b>
</pre></blockquote>

<p>
Interpretations of matrix <b>K</b>:
</p>

<ul>
<li> <b>K</b> is the value of the transfer function G(s) at s=0</li>
<li> For a stable state space system, a step input u[j] results in
     the output y[i](t->t<sub>&infin;</sub>) = K[i,j].</li>
</ul>

<p>
If <b>A</b> is singular (e.g. due to a zero eigenvalue), then a unique inverse
of <b>A</b> does not exist. If there are non-unique solutions of the
equation \"<b>A</b>*<b>X</b>=<b>B</b>\", the one with the smallest norm
in <b>X</b> is used to compute <b>K</b>. If no solution of this equation exists,
<b>K</b> cannot be computed. In this case, output argument
<b>finite</b> = <b>false</b> and all elements of
<b>K</b> are set to Modelica.Constants.inf (when <b>K</b> could be computed,
<b>finite</b> = <b>true</b>).
</p>
</html>"));
    end dcGain;

    encapsulated function isControllable
      "Check controllability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=
        Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

      output Boolean controllable "True, if ss is controllable";

    algorithm

      controllable := if StateSpace.Internal.isSISO(ss) then
        StateSpace.Internal.isControllableSISO(ss) else
        StateSpace.Internal.isControllableMIMO(ss, method);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
controllable = StateSpace.Analysis.<b>isControllable</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function StateSpace.Analysis.<b>isControllable</b> checks the controllability of a state space system. Therefore, the system is transformed into staircase form, i.e. the system matrix <b>H</b> of the transformed system has block upper Hessenberg form:
</p>
<blockquote><pre>
     | H11    H12     H13    ...     H1k |
     | H21    H22     H23    ...     H2k |
 <b>H</b> = |  0     H32     ...    ...     ... |
     | ...    ...     ...    ...     ... |
     |  0     ...      0    Hk,k-1   Hkk |
</pre></blockquote>
<p>
where, if <b>H</b>k,k-1 has full rank, indicating whether the system is controllable or not.
</p>
<p>
For single input systems the staircase form is a usual upper Hessenberg form, i.e. th blocks are of dimension one.<br>
The boolean input <b>method</b> defines for multi output systems the method to generate the staircase form of the system, whereas Types.StaircaseMethod.QR and Types.StaircaseMethod.SVD denotes QR-factorization and singular value decomposition respectively. Since staircase algorithm contains rank decisions QR-factorization should be restricted to well conditioned systems of lower order (&lt;5). Default is SVD.
</p>
<p>
Since controllability is dual to observability of the dual system (A', C', B', D'), proof of <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isObservable\">observability</a> is referred to proof of controllability of the dual system.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, -1, 1;0,0,-1],
    B=[0;0;1],
    C=[0,1,0],
    D=[0]);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean controllable;

<b>algorithm</b>
  controllable := Modelica_LinearSystems2.StateSpace.Analysis.isControllable(ss, method);
// controllable = true
</pre></blockquote>
</html>", revisions="<html>
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
    end isControllable;

    encapsulated function isObservable
      "Check observability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=
        Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

      output Boolean observable "True, if ss is observable";
    algorithm

      observable := if StateSpace.Internal.isSISO(ss) then
        StateSpace.Internal.isObservableSISO(ss) else
        StateSpace.Internal.isObservableMIMO(ss, method);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
observable = StateSpace.Analysis.<b>isObservable</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Function StateSpace.Analysis.<b>isObservable</b> checks the observability of a state space system. Since observability is dual to controllability of the dual system (A', C', B', D'), proof of observability is referred to proof of <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.isControllable\">controllability</a> of the dual system.<br>
The boolean input <b>method</b> defines for multi output systems the method to generate the staircase form of the system, whereas Types.StaircaseMethod.QR and Types.StaircaseMethod.SVD denotes QR-factorization and singular value decomposition respectively. Since staircase algorithm contains rank decisions QR-factorization should be restricted to  well conditioned systems of lower order (&lt;5). Default is SVD.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, -1, 1;0,0,-1],
    B=[0;0;1],
    C=[0,1,0],
    D=[0]);

  Types.Method method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD

  Boolean observable;

<b>algorithm</b>
  observable := Modelica_LinearSystems2.StateSpace.Analysis.isObservable(ss, method);
// observable = false
</pre></blockquote>
</html>", revisions="<html>
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
    end isObservable;

    encapsulated function isStabilizable
      "Check stabilizability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean stabilizable "True, if ss is stabilizable";

    algorithm
      if StateSpace.Internal.isSISO(ss) then
        stabilizable := StateSpace.Internal.isStabilizableSISO(ss);
      else
        stabilizable := StateSpace.Internal.isStabilizableMIMO(ss);
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
stabilizable = StateSpace.Analysis.<b>isStabilizable</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function checks whether a state space system is stabilizable or not.<br>
A system is stabilizable for the continuous-time case if all of the uncontrollable eigenvalues have negative real part.
Therefore, a controllable system is always stabilizable.
</p>
<p>
To check stabilizability, staircase algorithm is used to separate the controllable subspace from the uncontrollable subspace.
Then, the uncontrollable poles are checked to be stable, i.e. to have negative real parts.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Boolean stabilizable;

<b>algorithm</b>
   stabilizable := Modelica_LinearSystems2.StateSpace.Analysis.isStabilizable(ss);
// stabilizable = true
</pre></blockquote>
</html>", revisions="<html>
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
    end isStabilizable;

    encapsulated function isDetectable
      "Check detectability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean detectable "True, if ss is detectable";

    algorithm
      if StateSpace.Internal.isSISO(ss) then
        detectable := StateSpace.Internal.isDetectableSISO(ss);
      else
        detectable := StateSpace.Internal.isDetectableMIMO(ss);
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
detectable = StateSpace.Analysis.<b>isDetectable</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function checks whether a state space system is detectable or not.<br>
A system is detectable for the continuous-time case if all of the unobservable eigenvalues have negative real part.
Therefore, a observable system is always detectable.
</p>
<p>
To check detectability, staircase algorithm is used to separate the observable subspace from the unobservable subspace.
Then, the unobservable poles are checked to be stable, i.e. to have negative real parts.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, 1, 1;0,0,1],
    B=[0;0;1],
    C=[0,1,0],
    D=[0]);

  Boolean detectable;

<b>algorithm</b>
  detectable := Modelica_LinearSystems2.StateSpace.Analysis.isDetectable(ss);
// detectable = true
</pre></blockquote>
</html>", revisions="<html>
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
    end isDetectable;

    encapsulated function controllabilityMatrix
      "Calculate the controllability matrix [B, A*B, ..., A^(n-1)*B] of a state space system"

      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      output Real cm[size(ss.B, 1), size(ss.A, 2)*size(ss.B, 2)];

    algorithm
      if size(ss.A, 2) == 0 then
        cm := fill(
              0,
              size(ss.B, 1),
              0);
      else
        cm[:, 1:size(ss.B, 2)] := ss.B;

        for i in 2:size(ss.A, 1) loop
          cm[:, ((i - 1)*size(ss.B, 2) + 1):(i*size(ss.B, 2))] := ss.A*cm[:, ((
            i - 2)*size(ss.B, 2) + 1):((i - 1)*size(ss.B, 2))];
        end for;
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = StateSpace.Analysis.<b>controllabilityMatrix</b>(ss, method)
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

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Real Q[3,3];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.StateSpace.Analysis.controllabilityMatrix(ss);
// Q = [0, 1, 3; 0, 1, 2; 1, 1, 1]
</pre></blockquote>
</html>"));
    end controllabilityMatrix;

    encapsulated function observabilityMatrix
      "Calculate the observability matrix of a state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      output Real om[size(ss.A, 1)*size(ss.C, 1), size(ss.C, 2)];

    algorithm
      if size(ss.A, 1) == 0 then
        om := fill(
              0,
              0,
              size(ss.C, 2));
      else
        om[1:size(ss.C, 1), :] := ss.C;

        for i in 2:size(ss.A, 1) loop
          om[((i - 1)*size(ss.C, 1) + 1):(i*size(ss.C, 1)), :] := om[((i - 2)*
            size(ss.C, 1) + 1):((i - 1)*size(ss.C, 1)), :]*ss.A;
        end for;
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Q = StateSpace.Analysis.<b>observabilityMatrix</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the observability matrix
</p>
<blockquote>
  <b>Q</b> = [<b>C</b>; <b>C</b>*<b>A</b>; ...; <b>C</b>*<b>A</b>^(n-1)]
</blockquote>
<p>
of the system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
    <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Real Q[3,3];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.StateSpace.Analysis.observabilityMatrix(ss);
// Q = [0, 1, 0; 0, 1, 1; 1, 1, 2]
</pre></blockquote>
</html>"));
    end observabilityMatrix;

    function analysis2
      "Perform a system analysis based on the poles and zeros of the system"

      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Eigenvalue;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.Utilities.Plot;
      import DymolaCommands;

      input StateSpace ss;

      input Internal.AnalyseOptions analyseOptions=
          Modelica_LinearSystems2.Internal.AnalyseOptions(
              plotEigenValues=true,
              plotInvariantZeros=true,
              plotStepResponse=true,
              plotFrequencyResponse=true,
              printSystem=true,
              printEigenValues=true,
              printEigenValueProperties=true,
              printInvariantZeros=true,
              printControllability=true,
              printObservability=true,
              headingEigenValues="Eigenvalues",
              headingInvariantzeros="Invariant zeros",
              headingStepResponse="Step response",
              headingFrequencyResponse="Frequency response",
              dB_w = false);
      input String fileName="systemReport.html"
        "Name of html-file that contains eigenvalue table";
      input String systemName=""
        "Name of system (used as heading in html file)";
      input String description="" "Description of system (used in html file)";
    protected
      String dummyFileName="dummy" + fileName;
    public
      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

    protected
      StateSpace ssBalanced = StateSpace.Transformation.toBalancedForm(ss);
      Complex j = Modelica.ComplexMath.j;
      Eigenvalue ev[size(ss.A, 1)];
      Integer nx=size(ss.A, 1);
      Integer window=0;
      Real eval[nx, 2];
      Real revec[nx, nx];
      Real levec[nx, nx];
      Complex cev[size(ss.A, 1)];
      Complex systemZeros[:]=
          Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ssBalanced);
      Boolean isStable;
      Boolean isControllable;
      Boolean isStabilizable;
      Boolean isObservable;
      Boolean isDetectable;

      Real abs_evec[nx];
      String xNames2[nx];
      String heading="Eigenvalues" "Eigen values of system";
      Eigenvalue evSorted[size(ss.A, 1)];
      Integer evIndex[size(ss.A, 1)];
      Complex zerosSorted[:];
      Integer zerosIndex[size(systemZeros, 1)];
      Integer nReal;

      Integer i;
      Integer k;
      Complex evecComplex[size(ss.A, 1), size(ss.A, 1)];
      Plot.Records.Curve curves[2];
      Plot.Records.Diagram diagram2;
      Boolean instableZeros=false;

      String filePathOnly "NOT USED: Path to fileName";
      String fileNameOnly "Name of fileName without extension and path";
      String fileExtOnly "Extension of fileName";
      String fileNameImg "General name (without extension) of file for a plot";
      String fileNameImg2="none" "Current name of file for a plot";

      Internal.AnalyseOptions analyseOptions2=analyseOptions;
    algorithm
      // ---------------------------------------------------------------------------------------------------
      // The correct HTML format generated with this function can be checked with following commands:
      //   Modelica_LinearSystems2.StateSpace.Analysis.analysis(Modelica_LinearSystems2.StateSpace(A=[2,1,1;1,1,1;1,2,2], B=[1;2.2;3], C=[2,4,6;3,8,5], D=[6;4], yNames={"y1_test","y2_test"}, xNames={"xx1","xx2","xx3"}, uNames={"u1_test"}), description="Test file in HTML format from function 'analysis'.");
      //   Modelica_LinearSystems2.StateSpace.Analysis.analysis(Modelica_LinearSystems2.StateSpace(A=[2, 1.43, 12, 3; 1, 1, 1, 43; 1, 3, 2, 2; 1, 1, 4.2, 1.2], B=[1, 2; 2.2, 3; 3, 1; 4, 0], C=[25, 1.4, 6.3, 1; 0.3, 8, 5, 1; 1, 3, 2, 2], D=[6, 4; 4, 2; 6, 5], yNames={"y1_test","y2_te","y3_"}, xNames={"xx1","x2","xxx3","xx4"}, uNames={"u1_test","u2_test"}));
      // ---------------------------------------------------------------------------------------------------

      (filePathOnly,fileNameOnly,fileExtOnly) :=
        Modelica.Utilities.Files.splitPathName(fileName);
      fileNameImg := fileNameOnly;

      // If system has no inputs and outputs, modify analyze options that do not make sense
      if size(ss.B, 2) == 0 or size(ss.C, 1) == 0 then
        analyseOptions2.plotStepResponse := false;
        analyseOptions2.plotFrequencyResponse := false;
        analyseOptions2.printControllability := false;
        analyseOptions2.printObservability := false;
      end if;

      // If system has no states, modify analyze options that do not make sense
      if nx < 1 then
         analyseOptions2.plotEigenValues          :=false;
         analyseOptions2.plotInvariantZeros       :=false;
         analyseOptions2.printEigenValues         :=false;
         analyseOptions2.printEigenValueProperties:=false;
         analyseOptions2.printInvariantZeros      :=false;
      end if;

      // If system is too large, do not print A,B,C,D matrices
      if nx > 50 or size(ss.B, 2) > 50 or size(ss.C, 1) > 50 then
         analyseOptions2.printSystem:=false;
      end if;

      // Get eigenvalues
      // ---------------
      (eval,levec,revec) := Modelica_LinearSystems2.Math.Matrices.eigenValues(
        ss.A);

      for i in 1:nx loop
        cev[i].re := eval[i, 1];
        cev[i].im := eval[i, 2];
        ev[i].ev := cev[i];
      end for;

      (evSorted,evIndex) := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

      // Build x names
      // -------------
      if size(ss.xNames, 1) <> nx or nx > 1 and ss.xNames[1]=="" then
        for i in 1:nx loop
          xNames2[i] := "x[" + String(i) + "]";
        end for;
      else
        xNames2 := ss.xNames;
      end if;

      // Whole system checks
      // ===================
      // Stability check
      isStable := true;
      for i in 1:nx loop
        isStable := isStable and ev[i].ev.re < 0;
      end for;

      if size(ss.B, 2) > 0 and size(ss.C, 1) > 0 then
        // Controllability check, stabilizability check
        isControllable := StateSpace.Analysis.isControllable(ssBalanced);
        isStabilizable := StateSpace.Analysis.isStabilizable(ssBalanced);

        // Observability check, detectability check
        isObservable := StateSpace.Analysis.isObservable(ssBalanced);
        isDetectable := StateSpace.Analysis.isDetectable(ssBalanced);
      else
        isControllable := false;
        isStabilizable := false;
        isObservable := false;
        isDetectable := false;
      end if;

      // Analysis of single eingenvalues
      ev := StateSpace.Internal.characterizeEigenvalue(ss, ev);

      // Sort eigen values according to smallest imaginary value and restore the original order
      evSorted := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

      // Analysis file
      // -------------
      Modelica.Utilities.Files.removeFile(fileName);
      Modelica.Utilities.Files.removeFile(dummyFileName);

      // Text should be printed into new file in HTML environment
      // --------------------------------------------------------
      StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
      StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);

      if analyseOptions2.printSystem then
        printSystem(
              ss,
              fileName,
              systemName,
              description);

        printSystem(
              ss,
              dummyFileName,
              systemName,
              description);
      end if;

      printHead1(
            ss,
            isStable,
            isControllable,
            isStabilizable,
            isObservable,
            isDetectable,
            fileName,
            analyseOptions=analyseOptions2);
      printHead1(
            ss,
            isStable,
            isControllable,
            isStabilizable,
            isObservable,
            isDetectable,
            dummyFileName,
            analyseOptions=analyseOptions2);
      StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);

      Modelica.Utilities.Streams.readFile(dummyFileName);

      // Plot step response
      // ------------------
      if analyseOptions2.plotStepResponse then
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print(
          "<html>\n<body>\n<p>\n<b>Step responses</b>\n</p>\n</body>\n</html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);
        StateSpace.Plot.step(ss=ssBalanced);
        fileNameImg2 := fileNameImg + "Step.png";
        DymolaCommands.Plot.ExportPlotAsImage(fileName=fileNameImg2, id=-1, includeInLog=false);
        print("<p>\n<img src=\"" + fileNameImg2 + "\">\n</p>", fileName);
      end if;

      // Plot Bode plots
      if analyseOptions2.plotFrequencyResponse then
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print("<html>\n<body>\n<p>\n<b>Bode plots</b>\n</p>\n</body>\n</html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);
        StateSpace.Plot.bodeMIMO(ss=ss, Hz=not analyseOptions.dB_w, dB=analyseOptions.dB_w);
        //     fileNameImg2 := fileNameImg + "BodeMIMO1.png";
        //     ExportPlotAsImage(
        //       fileName=fileNameImg2,
        //       id=-1,
        //       includeInLog=false);
      end if;

      // Calculate the number of real eigenvalues
      nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(cev);

      // Construct complex eigenvector matrix
      i := 1;
      while i <= nx loop
        if eval[i, 2] == 0.0 then

          for jj in 1:nx loop
            evecComplex[jj, i] := revec[jj, i] + 0*j;
          end for;
          i := i + 1;
        else
          for jj in 1:nx loop
            evecComplex[jj, i] := revec[jj, i] + revec[jj, i + 1]*j;
            evecComplex[jj, i + 1] := revec[jj, i] - revec[jj, i + 1]*j;
          end for;
          i := i + 2;
        end if;
      end while;

      if analyseOptions2.printEigenValues then
        printHead2a(
              fileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal > 0));
        if nReal > 0 then
          printTab1(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                fileName,
                analyseOptions=analyseOptions2);
        end if;

        Modelica.Utilities.Files.removeFile(dummyFileName);
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead2a(
              dummyFileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal > 0));
        if nReal > 0 then
          printTab1(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                dummyFileName,
                analyseOptions=analyseOptions2);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
        Modelica.Utilities.Streams.readFile(dummyFileName);

        printHead2b(
              fileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal < nx));
        if nReal < nx then
          printTab2(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                fileName,
                analyseOptions=analyseOptions2);
        end if;

        Modelica.Utilities.Files.removeFile(dummyFileName);
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead2b(
              dummyFileName,
              analyseOptions=analyseOptions2,
              printTable=(nReal < nx));
        if nReal < nx then
          printTab2(
                evSorted,
                evIndex,
                revec,
                levec,
                nReal,
                xNames2,
                dummyFileName,
                analyseOptions=analyseOptions2);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
        Modelica.Utilities.Streams.readFile(dummyFileName);

       // Plot eigenvalues and invariant zeros
       if analyseOptions2.plotEigenValues or
          analyseOptions2.plotInvariantZeros and size(systemZeros, 1) > 0 then
          i := 0;
          if analyseOptions2.plotEigenValues then
            i := i + 1;
            curves[i] := Plot.Records.Curve(
                  x=eval[:, 1],
                  y=eval[:, 2],
                  legend="poles",
                  autoLine=false,
                  linePattern=Plot.Types.LinePattern.None,
                  lineSymbol=Plot.Types.PointSymbol.Cross);
          end if;

          // Plot invariant zeros
          if size(systemZeros, 1) > 0 and analyseOptions2.plotInvariantZeros then
            i := i + 1;
            curves[i] := Plot.Records.Curve(
                  x=systemZeros[:].re,
                  y=systemZeros[:].im,
                  legend="zeros",
                  autoLine=false,
                  linePattern=Plot.Types.LinePattern.None,
                  lineSymbol=Plot.Types.PointSymbol.Circle);
          end if;

          diagram2 := defaultDiagram;
          diagram2.curve := curves[1:i];
          Plot.diagram(diagram2, device);
        end if;

        if analyseOptions2.printEigenValueProperties then
          printHead3(fileName);
          printTab3(
                evSorted,
                evecComplex,
                evIndex,
                cev,
                nReal,
                xNames2,
                fileName);

          Modelica.Utilities.Files.removeFile(dummyFileName);
          StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
          printHead3(dummyFileName);
          printTab3(
                evSorted,
                evecComplex,
                evIndex,
                cev,
                nReal,
                xNames2,
                dummyFileName);
          StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);
          Modelica.Utilities.Streams.readFile(dummyFileName);

        end if;
      end if;

      // ZEROS
      (zerosSorted,zerosIndex) :=Modelica.ComplexMath.Vectors.sort(systemZeros);
      nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(zerosSorted);

      if analyseOptions2.printInvariantZeros then
        printHead4(fileName, printTable=(size(systemZeros, 1) > 0));
        if size(systemZeros, 1) > 0 then
          Modelica_LinearSystems2.StateSpace.Analysis.analysis.printTab4(
                zerosSorted,
                zerosIndex,
                nReal,
                fileName);
        end if;

        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, true);
        printHead4(dummyFileName, printTable=(size(systemZeros, 1) > 0));
        if size(systemZeros, 1) > 0 then
          printTab4(
                zerosSorted,
                zerosIndex,
                nReal,
                dummyFileName);
        end if;
        k := 0;
        for i in 1:size(systemZeros, 1) loop
          if systemZeros[i].re > 0 then
            k := k + 1;
          end if;
        end for;
        if k > 0 then
          print("<p>\n<b>Note, that the system has " + String(k) +
            " zeros in the right complex half-plane.</b>\n</p>", fileName);
          print("<p>\n<b>Note, that the system has " + String(k) +
            " zeros in the right complex half-plane.</b>\n</p>", dummyFileName);
        end if;
        StateSpace.Analysis.analysis.printHTMLbasics(dummyFileName, false);

      end if;
      Modelica.Utilities.Streams.readFile(dummyFileName);
      Modelica.Utilities.Files.removeFile(dummyFileName);

      print("\n\nAnalysis results have been written to file \"" +
        Modelica.Utilities.Files.fullPathName(fileName) + "\"");

      // Last print of HTML environment
      // --------------------------------------------------------
      StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);

      // SUB FUNCTIONS

    public
      encapsulated function printSystem
        "Print the state space system in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;
        import Modelica_LinearSystems2;

        input StateSpace ss "State space system to analyze";
        input String fileName="systemAnalysis.html"
          "File on which the state space system is written in html format";
        input String systemName="State Space System"
          "Name of the state space system";
        input String description="" "Description of system (used in html file)";
        input String format=".6g" "Format of numbers (e.g. \"20.8e\")";

        input Boolean htmlBasics=false
          "True, if text should be printed within 'html' and 'body' environment, otherwise text printed into existing file fileName"
          annotation (Dialog(group="HTML format"));
        input Integer hSize(
          min=1,
          max=5) = 1
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));
      protected
        Integer nx=size(ss.A, 1);
        Integer nu=size(ss.B, 2);
        Integer ny=size(ss.C, 1);
        Integer c1=integer(ceil(nx/2) - 1);
        Integer c2=integer(ceil(ny/2) - 1);
        Integer dist=2;
        Boolean centered=true
          "True, if matrices columns should be centered, otherwise right aligned";

        String td_align=if centered then "  <td style=\"text-align:center\">"
             else "  <td style=\"text-align:right\">";

      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 4 then 4 else hSize;
        String heading="h" + String(hSizeOK);
        String heading2="h" + String(hSizeOK + 1);
        String heading3="h" + String(hSizeOK + 2);
        Boolean printIndices;

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(Modelica_LinearSystems2.StateSpace(A=[2,1,1;1,1,1;1,2,2], B=[1;2.2;3], C=[2,4,6;3,8,5], D=[6;4], yNames={"y1_test","y2_test"}, xNames={"xx1","xx2","xx3"}, uNames={"u1_test"}));
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(Modelica_LinearSystems2.StateSpace(A=[2, 1.43, 12, 3; 1, 1, 1, 43; 1, 3, 2, 2; 1, 1, 4.2, 1.2], B=[1, 2; 2.2, 3; 3, 1; 4, 0], C=[25, 1.4, 6.3, 1; 0.3, 8, 5, 1; 1, 3, 2, 2], D=[6, 4; 4, 2; 6, 5], yNames={"y1_test","y2_te","y3_"}, xNames={"xx1","x2","xxx3","xx4"}, uNames={"u1_test","u2_test"}));
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printSystem(ss=Modelica_LinearSystems2.StateSpace.Import.fromModel("Modelica.Mechanics.Rotational.Examples.First"), description="Test file in HTML format from function printSystem.");
        // ---------------------------------------------------------------------------------------------------

        if htmlBasics then
          // Text should be printed into new file in HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
        end if;

        print("<" + heading + ">System report</" + heading + ">", fileName);
        print("\n<" + heading2 + ">General information</" + heading2 + ">",
          fileName);
        //print("<h1>System report</h1>", fileName);
        //print("\n<h2>General information</h2>", fileName);

        if systemName == "" then
        else
          print("\n<" + heading3 + ">System name</" + heading3 + ">", fileName);
          //print("\n<h3>System name</h3>", fileName);
          print("<p>\n" + systemName + "\n</p>", fileName);
        end if;

        if description == "" then
        else
          print("\n<" + heading3 + ">Description</" + heading3 + ">", fileName);
          //print("\n<h3>Description</h3>", fileName);
          print("<p>\n" + description + "\n</p>", fileName);
        end if;

        print("\n<" + heading3 + ">Matrices</" + heading3 + ">", fileName);
        //print("\n<h3>Matrices</h3>", fileName);
        print(
          "<p>\nThe system described in the state space representation\n</p>",
          fileName);
        print(
          "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin:20px 0 20px 20px;\" "
           + "cellpadding=\"3\" border=\"0\"> ", fileName);
        print("<tr><td>der(x) </td> <td>=</td> <td> Ax</td> <td> +</td><td> Bu</td></tr>
         <tr><td> y </td>     <td>=</td> <td> Cx</td> <td> + </td><td>Du</td></tr>",
          fileName);
        print("</table>\n<p>\nis defined by\n</p>", fileName);

        // ===============================
        // Print signal names and matrices (print row and column indices if at least one matrix has more as 5 elements)
        // ===============================
        printIndices := size(ss.A, 1) > 5 or size(ss.B, 2) > 5 or size(ss.C, 1)
           > 5;
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.uNames,
                "uNames",
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.yNames,
                "yNames",
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Vectors.printStringVectorInHtml(
                ss.xNames,
                "xNames",
                fileName=fileName,
                printIndices=printIndices);

        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.A,
                "A",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.B,
                "B",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.C,
                "C",
                format=format,
                fileName=fileName,
                printIndices=printIndices);
        Modelica_LinearSystems2.Math.Matrices.printMatrixInHtml(
                ss.D,
                "D",
                format=format,
                fileName=fileName,
                printIndices=printIndices);

        if ny == 0 and nu == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has neither inputs nor outputs (and therefore matrices B, C, and D are empty matrices)!\n</p>",
            fileName);
        elseif ny == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has no outputs (and therefore matrices C and D are empty matrices)!\n</p>",
            fileName);
        elseif nu == 0 then
          print(
            "<p>\n<b>Note</b>, that the system has no inputs (and therefore matrices B and D are empty matrices)!\n</p>",
            fileName);
        end if;

        if htmlBasics then
          // Last print of HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);
        end if;

      end printSystem;

      encapsulated function printHead1
        "Print the heading of document for characteristics in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input StateSpace ss;
        // This could be deleted sinc not used. But for reasons of beackward compatibility it is still here.
        input Boolean isStable;
        input Boolean isControllable;
        input Boolean isStabilizable;
        input Boolean isObservable;
        input Boolean isDetectable;
        input String fileName="systemHead1.html"
          "File on which the information is written in html format";

        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

        input Boolean htmlBasics=false
          "True, if text should be printed within 'html' and 'body' environment, otherwise text printed into existing file fileName"
          annotation (Dialog(group="HTML format"));
        input Integer hSize(
          min=1,
          max=5) = 2
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));

      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 5 then 5 else hSize;
        String heading="h" + String(hSizeOK);

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead1(Modelica_LinearSystems2.StateSpace(A=[2], B=[1], C=[1], D=[1]), false, false, false, true, false, htmlBasics=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        if htmlBasics then
          // Text should be printed into new file in HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, true);
        end if;

        print("\n<" + heading + ">Characteristics</" + heading +
          ">\n<p>\nThe system\n</p>\n<p> is ", fileName);

        if analyseOptions.printControllability and analyseOptions.printObservability then
          print((if isStable then " " else "<b>not</b> ") + "stable" + "\n<br>"
             + (if isStable then if isControllable then "and it is " else
            "but it is <b>not</b> " else if isControllable then "but it is "
             else "and it is <b>not</b> ") + "controllable" + (if isStable
             then "" else "\n<br>" + (if isControllable then
            " and therefore it is " else if isStabilizable then " but it is "
             else "and is <b>not</b> ") + "stabilizable.") +
            "\n<br> The system is " + (if isObservable then " " else
            "<b>not</b> ") + "observable" + (if isStable then "" else "\n<br>"
             + (if isObservable then " and therefore it is " else if
            isDetectable then " but it is " else "and is <b>not</b> ") +
            "detectable.") + "\n<br>", fileName);
        elseif not analyseOptions.printObservability and analyseOptions.printControllability then
          print((if isStable then " " else "<b>not</b> ") + "stable" + "\n<br>"
             + (if isStable then if isControllable then "and it is " else
            "but it is <b>not</b> " else if isControllable then "but it is "
             else "and it is <b>not</b> ") + "controllable" + (if isStable
             then "" else "\n<br>" + (if isControllable then
            " and therefore it is " else if isStabilizable then " but it is "
             else "and is <b>not</b> ") + "stabilizable.") + "\n<br>", fileName);
        elseif not analyseOptions.printControllability and analyseOptions.printObservability then
          print((if isStable then " " else "<b>not</b> ") + "stable." +
            "\n<br> The system is " + (if isObservable then " " else
            "<b>not</b> ") + "observable" + (if isStable then "" else "\n<br>"
             + (if isObservable then " and therefore it is " else if
            isDetectable then " but it is " else "and is <b>not</b> ") +
            "detectable.") + "\n<br>", fileName);
        else
          print((if isStable then " " else "<b>not</b> ") + "stable." +
            "\n<br>", fileName);
        end if;

        print("</p>", fileName);

        if htmlBasics then
          // Last print of HTML environment
          // --------------------------------------------------------
          StateSpace.Analysis.analysis.printHTMLbasics(fileName, false);
        end if;

      end printHead1;

      encapsulated function printHead2a
        "Print the heading of document for eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input String fileName="systemHead2a.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

        input Boolean printTable=true
          "True, if the system has real eigenvalues to be printed in table";
        input Integer hSize(
          min=1,
          max=5) = 3
          "Size of heading of printed document (=1: Title, =2: Chapter, etc.)"
          annotation (Dialog(group="HTML format"));
      protected
        Integer hSizeOK=if hSize < 1 then 1 else if hSize > 5 then 5 else hSize;
        String heading="h" + String(hSizeOK);

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead2a(htmlBasics=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        print("\n<" + heading + ">Eigenvalues analysis</" + heading + ">",
          fileName);
        //print("<p>\n<b>Real eigenvalues</b>\n</p>", fileName);

        if printTable then
          print("<p>\nThe system has the following real eigenvalues.\n</p>",
            fileName);
          print(
            "<table style=\"background-color:rgb(100, 100, 100);margin:20px 0 20px 20px;\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Real eigenvalues</caption>", fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> eigenvalue </td>\n  <td> T [s] </td>\n  <td> characteristics </td>",
            fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td> contribution to states</td>", fileName);
          end if;

          print("</tr>", fileName);
        else
          print("<p>\nThe system has no real eigenvalues.\n</p>", fileName);
        end if;

      end printHead2a;

      encapsulated function printHead2b
        "Print the heading of document for conjugated complex pairs in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

        input String fileName="systemHead2b.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");
        input Boolean printTable=true
          "True, if the system has complex pairs to be printed in table";

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead2b();
        // ---------------------------------------------------------------------------------------------------

        if printTable then
          print(
            "<p>\nThe system has the following complex conjugate pairs of eigenvalues.\n</p>",
            fileName);
          print(
            "<table style=\"background-color:rgb(100, 100, 100);margin:20px 0 20px 20px;\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Complex conjugate pairs of eigenvalues</caption>",
            fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> eigenvalue </td>\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n  <td> characteristics </td>",
            fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td> contribution to states</td>", fileName);
          end if;

          print("</tr>", fileName);
        else
          print(
            "<p>\nThe system has no complex conjugate eigenvalue pairs.\n</p>",
            fileName);
        end if;

      end printHead2b;

      encapsulated function printHead3
        "Print the heading of document for description in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

        input String fileName="systemHead3.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      algorithm
        print(
          "<p>\nIn the tables above, the column <b>contribution to states</b> lists for each eigenvalue the states to which the"
           +
          " corresponding modal state z[i] contributes most. This information is based on the"
           +
          " two largest absolute values of the corresponding right eigenvector (if the second large value"
           +
          " is less than 5&nbsp;% of the largest contribution, it is not shown). Note"
           +
          " the <b>right eigenvector</b> v<sub>j</sub> and the <b>left eigenvector</b> u<sub>j</sub> of A satisfy the"
           +
          " following relationships with regards to <b>eigenvalue</b> &lambda;<sub>j</sub>,"
           +
          " state vector x and modal state vector z (u<sub>j</sub><sup>H</sup> denotes the conjugate transpose of u<sub>j</sub>):"
           + " </p>" +
          " <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\">" +
          " <tr><td width=\"50\"></td>" +
          "\n    <td>A * v<sub>j</sub> = &lambda;<sub>j</sub> * v<sub>j</sub>; &nbsp;&nbsp;&nbsp;&nbsp;"
           +
          "         u<sub>j</sub><sup>H</sup> * A = &lambda;<sub>j</sub> * u<sub>j</sub><sup>H</sup>; &nbsp;&nbsp;&nbsp;&nbsp;"
           +
          "               x = V * z; &nbsp;&nbsp;&nbsp;&nbsp; V = [v<sub>1</sub>, v<sub>2</sub>, ...]</td>"
           + "           </tr>" + "\n</table>" + "\n<p>" +
          "\nIn the next table, for each state in the column <b>correlation to modal states</b>, the modal"
           +
          " states z[i] which contribute most to the corresponding state are summarized, that is"
           + " the state is mostly composed of these modal states." +
          "\nThis information is based on the two largest absolute values of row i of the"
           +
          " eigenvector matrix that is associated with eigenvalue i (if the second large value"
           +
          " is less than 5&nbsp;% of the largest contribution, it is not shown). This only holds"
           +
          " if the modal states z[i] are in the same order of magnitude. Otherwise, the listed modal states"
           + " might be not the most relevant ones.</p>", fileName);

        print(
          "<table style=\"background-color:rgb(100, 100, 100); margin:20px 0 20px 20px;\" "
           + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">\n" +
          "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           +
          "\n  <td> state </td>\n  <td> correlation to modal states </td>\n  <td> eigenvalue # </td>"
           +
          "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n  <td> T [s] </td>\n</tr>",
          fileName);

      end printHead3;

      encapsulated function printHead4
        "Print the heading of document for invariant zeros in html format on file"
        import Modelica;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

        input String fileName="systemHead4.html"
          "File on which the information is written in html format";
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");
        input Boolean printTable=true
          "True, if the system has complex pairs to be printed in table";

      algorithm
        // ---------------------------------------------------------------------------------------------------
        // The correct HTML format generated with this function can be checked with following commands:
        //   Modelica_LinearSystems2.StateSpace.Analysis.analysis.printHead4(htmlEnv=true, hSize=3);
        // ---------------------------------------------------------------------------------------------------

        if printTable then
          print("<p>\nThe system has the following invariant zeros.\n</p>",
            fileName);
          print(
            "\n<table style=\"background-color:rgb(100, 100, 100); margin:20px 0 20px 20px;\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", fileName);
          print("<caption>Invariant zeros</caption>", fileName);
          print(
            "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
             +
            "\n  <td> number </td>\n  <td> invariant zero </td>\n  <td> Time constant [s] </td>"
             + "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n</tr>",
            fileName);
        else
          print("<p>\nThe system has no invariant zeros.\n</p>", fileName);
        end if;

      end printHead4;

      encapsulated function printHTMLbasics
        "Print the html preamble or ending on file"
        import Modelica.Utilities.Files;
        import Modelica.Utilities.Streams;

        input String fileName="systemReport.html"
          "File on which the html basics should be written";
        input Boolean printBegin=false
          "True, if beginning of a html file should be printed, otherwise the ending"
          annotation (choices(checkBox=true));

      algorithm
        if printBegin then
          // First print of HTML environment into new file
          Files.removeFile(fileName);
          // Following doesn't work in Dymola
          //Streams.print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", fileName);
          Streams.print("<html>", fileName);
          Streams.print(
            "<head>\n  <title>Analysis of a state space system from Modelica LinearSystems2</title>\n</head>",
            fileName);
          Streams.print("<style type=\"text/css\">", fileName);
          Streams.print("* { font-size: 10pt; font-family: Arial,sans-serif; }",
            fileName);
          Streams.print("</style>", fileName);
        else
          // Last print of HTML environment
          Streams.print("</html>", fileName);
        end if;
      end printHTMLbasics;

      encapsulated function printTab1
        "Print the table with real eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Complex;

        input Eigenvalue evSorted[:];
        input Integer evIndex[size(evSorted, 1)];
        input Real r_evec[size(evSorted, 1), size(evSorted, 1)];
        input Real l_evec[size(evSorted, 1), size(evSorted, 1)];
        input Integer nReal;
        input String xNames2[size(evSorted, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evSorted, 1);
        Real w;
        Real d;

        Integer i;
        Integer j;
        Integer k;
        String number;

        Real r_abs_evec[nx];
        Real l_abs_evec[nx];
        Integer r_maxIndex1;
        Integer l_maxIndex1;
        Integer r_maxIndex2;
        Integer l_maxIndex2;
        //  Complex v_normalized[size(evSorted,1)];
        Real r_abs_v_normalized;
        Real l_abs_v_normalized;
        Real r_v;
        Real l_v;
        Real r_absMax1;
        Real l_absMax1;
        Real r_absMax2;
        Real l_absMax2;
        Boolean r_two;
        Boolean l_two;
        Boolean r_first;
        Boolean l_first;

      algorithm
        i := 1;
        j := i;
        while i <= nReal loop
          // Build eigenvalue number

          number := String(
                  i,
                  minimumLength=7,
                  leftJustified=false);
          j := j + 1;

          // Determine largest value in eigenvector
          k := evIndex[i] "Index with respect to unsorted eigen values";
          r_abs_evec := abs(r_evec[:, k]);
          l_abs_evec := abs(l_evec[:, k]);

          r_first := true;
          r_two := false;
          r_absMax1 := 0;
          r_maxIndex1 := 0;
          r_absMax2 := 0;
          r_maxIndex2 := 0;
          r_abs_v_normalized := Modelica.Math.Vectors.norm(r_abs_evec, 1);

          l_first := true;
          l_two := false;
          l_absMax1 := 0;
          l_maxIndex1 := 0;
          l_absMax2 := 0;
          l_maxIndex2 := 0;
          l_abs_v_normalized := Modelica.Math.Vectors.norm(l_abs_evec, 1);
          for j in 1:nx loop
            r_v := r_abs_evec[j];
            l_v := l_abs_evec[j];

            if r_first then
              r_first := false;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif not r_two then
              r_two := true;
              if r_v < r_absMax1 then
                r_absMax2 := r_v;
                r_maxIndex2 := j;
              else
                r_absMax2 := r_absMax1;
                r_maxIndex2 := r_maxIndex1;
                r_absMax1 := r_v;
                r_maxIndex1 := j;
              end if;
            elseif r_v > r_absMax1 then
              r_absMax2 := r_absMax1;
              r_maxIndex2 := r_maxIndex1;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif r_v > r_absMax2 then
              r_absMax2 := r_v;
              r_maxIndex2 := j;
            end if;

            if l_first then
              l_first := false;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif not l_two then
              l_two := true;
              if l_v < l_absMax1 then
                l_absMax2 := l_v;
                l_maxIndex2 := j;
              else
                l_absMax2 := l_absMax1;
                l_maxIndex2 := l_maxIndex1;
                l_absMax1 := l_v;
                l_maxIndex1 := j;
              end if;
            elseif l_v > l_absMax1 then
              l_absMax2 := l_absMax1;
              l_maxIndex2 := l_maxIndex1;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif l_v > l_absMax2 then
              l_absMax2 := l_v;
              l_maxIndex2 := j;
            end if;

          end for;

          r_absMax1 := 100*r_absMax1/r_abs_v_normalized;
          r_absMax2 := 100*r_absMax2/r_abs_v_normalized;

          if r_absMax2 < 0.05*r_absMax1 then
            r_two := false;
          end if;

          l_absMax1 := 100*l_absMax1/l_abs_v_normalized;
          l_absMax2 := 100*l_absMax2/l_abs_v_normalized;

          if l_absMax2 < 0.05*l_absMax1 then
            l_two := false;
          end if;

          // Print data for one eigen value
          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:center\"> "
             + number + " </td>\n  <td style=\"text-align:left\"> &nbsp; " +
            String(evSorted[i].ev.re, format="14.4e") +
            " </td>\n  <td style=\"text-align:left\"> &nbsp; " + (if evSorted[i].timeConstant
             < 1e6 then String(evSorted[i].timeConstant, format="9.4f") else
            "---") + " </td>\n  <td style=\"text-align:left\"> &nbsp; " + (if
            evSorted[i].isStable then "" else "not ") + "stable, " + (if
            evSorted[i].isStable then (if evSorted[i].isControllable then ""
             else "not ") + "controllable, " else (if evSorted[i].isStabilizable
             then "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable
             then (if evSorted[i].isObservable then "" else "not ") +
            "observable " else (if evSorted[i].isDetectable then "" else "not ")
             + "detectable ") + " </td>", fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td style=\"text-align:left\"> &nbsp; " + " z[" + String(i)
               + "]" + " contributes to " + xNames2[r_maxIndex1] + " with " +
              String(r_absMax1, format=".3g") + " %<br>" + (if r_two then
              "&nbsp; " + " z[" + String(i) + "]" + " contributes to " +
              xNames2[r_maxIndex2] + " with " + String(r_absMax2, format=".3g")
               + " %" else "") + " </td>", fileName);
          end if;

          print("</tr>", fileName);

          i := j;
        end while;

        print("</table>", fileName);
      end printTab1;

      encapsulated function printTab2
        "Print the table with complex conjugate eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Complex;

        input Eigenvalue evSorted[:];
        input Integer evIndex[size(evSorted, 1)];
        input Real r_evec[size(evSorted, 1), size(evSorted, 1)];
        input Real l_evec[size(evSorted, 1), size(evSorted, 1)];
        input Integer nReal;
        input String xNames2[size(evSorted, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evSorted, 1);

        Integer i;
        Integer k;
        String number;
        String number2;
        Integer j;

        Real r_abs_evec[nx];
        Real l_abs_evec[nx];
        Integer r_maxIndex1;
        Integer r_maxIndex2;
        Integer l_maxIndex1;
        Integer l_maxIndex2;
        Real r_abs_v_normalized;
        Real l_abs_v_normalized;
        Real r_v;
        Real l_v;
        Real r_absMax1;
        Real r_absMax2;
        Real l_absMax1;
        Real l_absMax2;
        Boolean r_two;
        Boolean l_two;
        Boolean r_first;
        Boolean l_first;

      algorithm
        i := nReal + 1;
        j := i;
        while i <= nx loop
          // Build eigenvalue number
          number := String(i) + "/" + String(i + 1);
          number2 := number;
          number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;
          j := j + 2;

          // Determine largest value in eigenvector
          k := evIndex[i] "Index with respect to unsorted eigen values";

          for i2 in 1:nx loop
            r_abs_evec[i2] := sqrt(r_evec[i2, k]^2 + r_evec[i2, k + 1]^2);
            l_abs_evec[i2] := sqrt(l_evec[i2, k]^2 + l_evec[i2, k + 1]^2);
          end for;

          r_first := true;
          r_two := false;
          r_absMax1 := 0;
          r_maxIndex1 := 0;
          r_absMax2 := 0;
          r_maxIndex2 := 0;
          r_abs_v_normalized := Modelica.Math.Vectors.norm(r_abs_evec, 1);
          l_first := true;
          l_two := false;
          l_absMax1 := 0;
          l_maxIndex1 := 0;
          l_absMax2 := 0;
          l_maxIndex2 := 0;
          l_abs_v_normalized := Modelica.Math.Vectors.norm(l_abs_evec, 1);

          for j in 1:nx loop
            r_v := r_abs_evec[j];
            l_v := l_abs_evec[j];

            if r_first then
              r_first := false;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif not r_two then
              r_two := true;
              if r_v < r_absMax1 then
                r_absMax2 := r_v;
                r_maxIndex2 := j;
              else
                r_absMax2 := r_absMax1;
                r_maxIndex2 := r_maxIndex1;
                r_absMax1 := r_v;
                r_maxIndex1 := j;
              end if;
            elseif r_v > r_absMax1 then
              r_absMax2 := r_absMax1;
              r_maxIndex2 := r_maxIndex1;
              r_absMax1 := r_v;
              r_maxIndex1 := j;
            elseif r_v > r_absMax2 then
              r_absMax2 := r_v;
              r_maxIndex2 := j;
            end if;

            if l_first then
              l_first := false;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif not l_two then
              l_two := true;
              if l_v < l_absMax1 then
                l_absMax2 := l_v;
                l_maxIndex2 := j;
              else
                l_absMax2 := l_absMax1;
                l_maxIndex2 := l_maxIndex1;
                l_absMax1 := l_v;
                l_maxIndex1 := j;
              end if;
            elseif l_v > l_absMax1 then
              l_absMax2 := l_absMax1;
              l_maxIndex2 := l_maxIndex1;
              l_absMax1 := l_v;
              l_maxIndex1 := j;
            elseif l_v > l_absMax2 then
              l_absMax2 := l_v;
              l_maxIndex2 := j;
            end if;

          end for;
          r_absMax1 := 100*r_absMax1/r_abs_v_normalized;
          r_absMax2 := 100*r_absMax2/r_abs_v_normalized;
          if r_absMax2 < 0.05*r_absMax1 then
            r_two := false;
          end if;

          l_absMax1 := 100*l_absMax1/l_abs_v_normalized;
          l_absMax2 := 100*l_absMax2/l_abs_v_normalized;
          if l_absMax2 < 0.05*l_absMax1 then
            l_two := false;
          end if;

          // Print data for one eigen value
          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> "
             + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " &plusmn; " +
            String(evSorted[i].ev.im, format="12.4e") + "j" + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(evSorted[i].frequency,
            format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(evSorted[i].damping,
            format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + (if evSorted[i].isStable
             then "" else "not ") + "stable, " + (if evSorted[i].isStable then
            (if evSorted[i].isControllable then "" else "not ") +
            "controllable, " else (if evSorted[i].isStabilizable then "" else
            "not ") + "stabilizable, ") + (if evSorted[i].isStable then (if
            evSorted[i].isObservable then "" else "not ") + "observable " else
            (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td>", fileName);

          if analyseOptions.printEigenValueProperties then
            print("  <td style=\"text-align:left\"> &nbsp; " + " z[" + number2
               + "]" + " contribute to " + xNames2[r_maxIndex1] + " with " +
              String(r_absMax1, format=".3g") + " %<br>" + (if r_two then
              "&nbsp; " + " z[" + number2 + "]" + " contribute to " + xNames2[
              r_maxIndex2] + " with " + String(r_absMax2, format=".3g") + " %"
               else "") + " </td>", fileName);
          end if;

          print("</tr>", fileName);

          i := j;
        end while;

        print("</table>", fileName);
      end printTab2;

      encapsulated function printTab3
        "Print the table with eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Complex;

        input Eigenvalue evSorted[:];
        input Complex evecComplex[:, :];
        input Integer evIndex[size(evecComplex, 1)];
        input Complex cev[size(evecComplex, 1)];
        input Integer nReal;
        input String xNames2[size(evecComplex, 1)];
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nx=size(evecComplex, 1);

        Integer maxIndex1;
        Integer maxIndex2;

        Complex v_normalized[size(evecComplex, 1)];
        Real abs_v_normalized;
        Real v;
        Real absMax1;
        Real absMax2;
        Boolean two;
        Boolean first;
        Integer j;
        Integer k;
        Integer iw1;
        Integer iw2;
        String number1;
        String number2;
        Real w1;
        Real w2;
        Real d1;
        Real d2;

      algorithm
        for i in 1:nx loop
          // Normalize i-th row of complex eigenvector matrix and determine two largest elements
          v_normalized :=Modelica.ComplexMath.Vectors.normalize(evecComplex[i, :]);
          first := true;
          two := false;
          absMax1 := 0;
          maxIndex1 := 0;
          absMax2 := 0;
          maxIndex2 := 0;
          j := 1;
          abs_v_normalized :=Modelica.ComplexMath.Vectors.norm(v_normalized, 1);
          while j <= nx loop
            if cev[j].im == 0 then
              v := abs(v_normalized[j].re);
              k := j;
              j := j + 1;
            else
              v :=2*Modelica.ComplexMath.'abs'(v_normalized[j]);
              k := j;
              j := j + 2;
            end if;

            if first then
              first := false;
              absMax1 := v;
              maxIndex1 := k;
            elseif not two then
              two := true;
              if v < absMax1 then
                absMax2 := v;
                maxIndex2 := k;
              else
                absMax2 := absMax1;
                maxIndex2 := maxIndex1;
                absMax1 := v;
                maxIndex1 := k;
              end if;
            elseif v > absMax1 then
              absMax2 := absMax1;
              maxIndex2 := maxIndex1;
              absMax1 := v;
              maxIndex1 := k;
            elseif v > absMax2 then
              absMax2 := v;
              maxIndex2 := k;
            end if;
          end while;

          if abs_v_normalized > 1e-30 then
            absMax1 := absMax1/abs_v_normalized;
            absMax2 := absMax2/abs_v_normalized;
          end if;

          if absMax2 < 0.05*absMax1 then
            two := false;
          end if;

          // Determine frequency and number of corresponding eigenvalue
          (w1,d1) :=Modelica_LinearSystems2.Math.ComplexAdvanced.frequency(cev[maxIndex1]);
          iw1 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex1, evIndex);
          if iw1 <= nReal then
            number1 := String(iw1);
          else
            number1 := String(iw1) + "/" + String(iw1 + 1);
          end if;

          if two then
            (w2,d2) :=Modelica_LinearSystems2.Math.ComplexAdvanced.frequency(cev[maxIndex2]);
            iw2 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex2, evIndex);
            if iw2 <= nReal then
              number2 := String(iw2);
            else
              number2 := String(iw2) + "/" + String(iw2 + 1);
            end if;
          end if;

          if two then
            print(
              "<tr style=\"background-color:white\">\n  <td rowspan=2 style=\"text-align:left\"> &nbsp; "
               + xNames2[i] + " </td>" +
              "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax1, format="5.1f") + "% by z[" + number1 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number1 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(w1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(d1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then String(evSorted[i].timeConstant, format="9.4f") else
              "--- </td>") + "\n</tr>\n<tr style=\"background-color:white\">"
               + "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax2, format="5.1f") + "% by z[" + number2 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number2 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw2 <= nReal
               then "---" else String(w2, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw2 <= nReal
               then "---" else String(d2, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if (iw2 <=
              nReal and abs(cev[maxIndex2].re) > 1e-10) then String(1/abs(cev[
              maxIndex2].re), format="9.4f") else "--- </td>\n</tr>"), fileName);
          else
            print(
              "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
               + xNames2[i] + " </td>" +
              "\n  <td style=\"text-align:left\"> &nbsp; is composed of " +
              String(100*absMax1, format="5.1f") + "% by z[" + number1 +
              "]</td>" + "\n  <td style=\"text-align:center\"> &nbsp; " +
              number1 + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(w1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then "---" else String(d1, format="9.4f")) + "</td>" +
              "\n  <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal
               then String(evSorted[i].timeConstant, format="9.4f") else
              "--- </td>\n</tr>"), fileName);
          end if;
          //     print("<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; " + xNames2[i] + " </td>\n  <td style=\"text-align:left\"> &nbsp; "
          //        + " is composed of " + String(100*absMax1, format="5.1f") + "% by z[" +
          //       number1 + "]" + (if two then " <br>" + " &nbsp; " + " is composed of " +
          //       String(100*absMax2, format="5.1f") + "% by z[" + number2 + "]" else "") + " </td> <td style=\"text-align:center\"> &nbsp; "
          //        + number1 + (if two then "<br> &nbsp; " + number2 else Strings.repeat(9))
          //        + " </td> <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal then
          //             "---" else String(w1, format="9.4f")) + (if two then "<br> &nbsp; "
          //        + (if iw2 <= nReal then "---" else String(w2, format="9.4f")) else
          //       Strings.repeat(9)) + " </td>\n  <td style=\"text-align:center\"> &nbsp; " +
          //       (if iw1 <= nReal then "---" else String(d1, format="9.4f")) + (if two then
          //             "<br> &nbsp; " + (if iw2 <= nReal then "---" else String(d2,
          //       format="9.4f")) else "") + " </td>\n  <td style=\"text-align:center\"> &nbsp; "
          //        + (if (iw1 <= nReal) then String(evSorted[i].timeConstant, format="9.4f") else
          //             "---") + (if two then "<br> &nbsp; " + (if (iw2 <= nReal and abs(
          //       cev[maxIndex2].re) > 1e-10) then String(1/abs(cev[maxIndex2].re),
          //       format="9.4f") else "---") else "") + " </td>\n</tr> ", fileName);

        end for;
        print("</table>", fileName);

      end printTab3;

      encapsulated function printTab4
        "Print the table with eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica.Utilities.Streams.print;
        import Complex;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;

        input Complex systemZeros[:];
        input Integer evIndex[size(systemZeros, 1)];
        input Integer nReal;
        input String fileName;
        input Modelica_LinearSystems2.Internal.AnalyseOptions analyseOptions=
            Modelica_LinearSystems2.Internal.AnalyseOptions(
                  plotEigenValues=true,
                  plotInvariantZeros=true,
                  plotStepResponse=true,
                  plotFrequencyResponse=true,
                  printEigenValues=true,
                  printEigenValueProperties=true,
                  printInvariantZeros=true,
                  printControllability=true,
                  printObservability=true,
                  headingEigenValues="Eigenvalues",
                  headingInvariantzeros="Invariant zeros",
                  headingStepResponse="Step response",
                  headingFrequencyResponse="Frequency response");

      protected
        Integer nz=size(systemZeros, 1);

        String number;
        Real timeConstant;
        Real freq;
        Real damp;

      algorithm
        for i in 1:nReal loop
          // Build eigenvalue number

          number := String(
                  i,
                  minimumLength=7,
                  leftJustified=false);
          timeConstant := if abs(systemZeros[i].re) > 10*Modelica.Constants.eps
             then 1/abs(systemZeros[i].re) else 1/(10*Modelica.Constants.eps);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td> &nbsp; " + String(systemZeros[i].re,
            format="14.4e") + " </td>" + "\n  <td> &nbsp; " + String(
            timeConstant, format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>\n</tr>",
            fileName);

        end for;

        for i in nReal + 1:2:nz loop
          number := String(i) + "/" + String(i + 1);
          number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;

          // Determine frequency and number of corresponding zero
          (freq,damp) :=Modelica_LinearSystems2.Math.ComplexAdvanced.frequency(systemZeros[i]);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(systemZeros[i].re, format="14.4e") + " &plusmn; " +
            String(systemZeros[i].im, format="12.4e") + "j </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(freq, format=
            "9.4f") + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(damp, format="9.4f") + " </td>\n</tr>", fileName);

        end for;

        print("</table>\n", fileName);
      end printTab4;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(ss);
   or
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(
  ss,
  analyseOptions=<a href=\"modelica://Modelica_LinearSystems2.Internal.AnalyseOptions\">analyseOptions</a>,
  fileName,
  systemName,
  description);
</pre></blockquote>

<h4>Description</h4>
<p>
This function analyzes a state space system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>     <label for=\"eqn1\">(1)</label>
    <b>x</b>(t=0) = <b>x</b><sub>0</sub>
</pre></blockquote>
<p>
based on its poles, i.e. the eigenvalues, and the zeros of the system.
The system will be checked for stability, controllability and observability. In the case that the system is not stable stabilizability and detectability are examined. Furthermore, stability, controllability, observability, stabilizability, and detectability are indicated for each eigenvalue.
</p>

<h5>Stability</h5>
<p>
System (1) is stable if and only if all eigenvalues of the matrix <b>A</b> have negative real parts.
The calculation of the eigenvalues is based on the LAPACK routine dgeev.
</p>

<h5>Controllability</h5>
<p>
System (1) is said to be controllable if, starting from any initial state <b>x</b><sub>0</sub>, the system can be driven by appropriate inputs to any final state <b>x</b><sub>1</sub> within some finite time window. Equivalent is that the eigenvalues of <b>A</b>-<b>BK</b> can  arbitrarily be assigned by an appropriate choice of the matrix <b>K</b>.
</p>

<h5>Stabilizability</h5>
<p>
System (1) is said to be stabilizable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are controllable. Therefore, a controllable system is always stabilizable. An equivalent definition of stabilizability is, that a system is said to be stabilizable if there exist a matrix <b>K</b> such that <b>A</b>-<b>BK</b> is stable.
</p>

<h5>Observability</h5>
<p>
System (1) is said to be observable if the (arbitrary) initial state <b>x</b><sub>0</sub> can be uniquely determined from any state <b>x</b>(t<sub>1</sub>), t<sub>1</sub>>0, from the knowledge of the input <b>u</b>(t) and output <b>y</b>(t). With other words,  from the system's outputs it is possible to determine the behavior of the entire system. Equivalent is, that the eigenvalues of <b>A</b>-<b>LC</b> can be arbitrarily be assigned by an appropriate choice of matrix <b>L</b>.
Observability is called the dual concept of controllability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is observable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is controllable.
</p>

<h5>Detectability</h5>
<p>
System (1) is said to be detectable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are observable. Therefore, a observable system is always detectable. An equivalent definition of detectability is, that a system is said to be detectable if there exist a matrix <b>L</b> such that <b>A</b>-<b>LC</b> is stable.
Detectability is called the dual concept of stabilizability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is detectable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is stabilizable.
</p>

<h5>Algorithm to test controllability/stabilizability and observability/detectability respectively</h5>
<p>
The test of controllability and stabilizability is performed with the staircase algorithm which transforms the system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) into the controller-Hessenberg form (<b>A</b><sub>H</sub>, <b>B</b><sub>H</sub>, <b>C</b><sub>H</sub>, <b>D</b>) with <b>A</b><sub>H</sub> is a block upper Hessenberg matrix and <b>B</b><sub>H</sub>=[<b>B</b><sub>1</sub>; 0] with triangular matrix <b>B</b><sub>1</sub> with rank(<b>B</b><sub>1</sub>) = rank(<b>B</b>).
In <b>A</b><sub>H</sub>=[<b>A</b><sub>c</sub>, *,0, <b>A</b><sub>nc</sub>) the eigenvalues of the matrices <b>A</b><sub>c</sub> and <b>A</b><sub>nc</sub> are the controllable eigenvalues and uncontrollable eigenvalues of <b>A</b> respectively.
The test of observability and detectability is performed by testing the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) with respect to controllability and stabilizability.
</p>

<h5>Solution of a linear time invariant system </h5>
<p>
The solution <b>x</b>(t) of the initial value problem (1) consists of the homogeneous part (zero input response) <b>x</b><sub>h</sub>(t) and the inhomogeneous part x<sub>i</sub>(t). The zero input solution is given by
</p>
<blockquote><pre>
<b>x</b><sub>h</sub>(t) = exp(<b>A</b>*(t-t<sub>0</sub>))<b>x</b><sub>0</sub>.
</pre></blockquote>
<p>
The system can also be represented as a linear combination of the modal states <b>z</b>,
</p>
<blockquote><pre>
<b>x</b> = <b>V</b><b>z</b>
</pre></blockquote>
<p>
i.e. the states of a similar system, with
</p>
<blockquote><pre>
der(<b>z</b>) = <b>V</b><sup>-1</sup><b>AVz</b> + <b>V</b><sup>-1</sup><b>B</b><b>u</b>
</pre></blockquote>
<p>
where the system matrix <b>V</b><sup>-1</sup><b>AV</b> is the real Jordan form. For single real eigenvectors the system is decoupled, i.e. the solution of the modal states are denoted by
<blockquote><pre>
z<sub>i</sub> = exp(s<sub>i</sub> t)*z<sub>0i</sub>
</pre></blockquote>
<p>
The behavior of the modal states is determined as the solution of a linear first order differential equation for real eigenvalues. Since this behavior is well known, the behavior of the x<sub>i</sub> can at least roughly be estimated by means of the behavior of the most relevant modal states. Therefore, the contribution of the modal states to the states is computed as an indication of the original system behavior.
</p>

<h5>Contribution of the modal states to the states</h5>
<p>
Generally, as described above, the states of the system can be described as linear combination of modal states and, therefore, the states can be characterized to a certain extend by the modal states if the proportions of the combination are known. Hence, for each modal state z<sub>i</sub> of the vector <b>z</b> the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub>| of the corresponding right eigenvector <b>v</b><sub>i</sub> indicate the proportion of <b>z</b><sub>i</sub> that is contributed to the state x<sub>j</sub>.
On the other hand, the composition of xi is indicated by the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>|, i.e. the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>| of the corresponding row <b>v</b><sub>i</sub><sup>T</sup> of the eigenvector matrix <b>V</b> indicate the proportion of the state x<sub>i</sub> that is contributed by the modal state z<sub>j</sub>.
</p>

<h4>Example</h4>
<blockquote><pre>
  ss=StateSpace(
    A=[-3,2,-3,4,5,6; 0,6,7,8,9,4; 0,2,3,0,78,6; 0,1,2,2,3,3; 0,13,34,0,0,1; 0,
      0,0,-17,0,0],
    B=[1,0; 0,1; 1,0; 0,1; 1,0; 0,1],
    C=[0,0,1,0,1,0; 0,1,0,0,1,1],
    D=[0,0; 0,0],
    xNames={\"x1\",\"x2\",\"x3\",\"x4\",\"x5\",\"x6\"},
    uNames={\"u1\",\"u2\"}, yNames={\"y1\",\"y2\"});

  String fileName=\"analysis.html\";
  String systemName=\"Demonstration System\";
  String description=\"System to demonstrate the usage of Modelica_LinearSystems2.StateSpace.Analysis.anlysis()\"

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Analysis.analysis(ss, fileName=fileName, systemName=systemName, description=description)
//  gives:
</pre></blockquote>

<h4>System report</h4>
<p>
The system <b>Demonstation System</b>
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>
</pre></blockquote>
<p>
is defined by
</p>
<blockquote><pre>
        x1   x2   x3   x4   x5   x6            u1  u2
    x1  -3    2   -3    4    5    6         x1  1   0
    x2   0    6    7    8    9    4         x2  0   1
A = x3   0    2    3    0   78    6     B = x3  1   0
    x4   0    1    2    2    3    3         x4  0   1
    x5   0   13   34    0    0    1         x5  1   0
    x6   0    0    0  -17    0    0         x6  0   1

        x1   x2   x3   x4   x5   x6            u1  u2
C = y1   0    0    1    0    1    0     D = y1  0   0
    y2   0    1    0    0    1    1         y2  0   0
</pre></blockquote>

<h5>Description</h5>
<p>
System to demonstrate the usage of Modelica_LinearSystems2.StateSpace.Analysis.analysis()
</p>

<h5>Characteristics</h5>
<p>The system
<br> is
not stable
<br>but it is controllable
<br> and therefore it is stabilizable
<br> The system is not observable
<br> but it is detectable
</p>

<p>
<b><big>Eigenvalues analysis</big></b>
<br><br>
<b>Real eigenvalues</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin:20px 0 20px 20px;\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td><td> eigenvalue </td> <td> T [s] </td>  <td> characteristics </td><td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:center\">       1 </td> <td style=\"text-align:left\"> &nbsp;   -4.9874e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0201 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, observable  </td> <td style=\"text-align:left\"> &nbsp;  z[1] contributes to x3 with 54.6 %<br>&nbsp;  z[1] contributes to x5 with 37 % </td> </tr>
<tr>
 <td style=\"text-align:center\">       2 </td> <td style=\"text-align:left\"> &nbsp;   -3.0000e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3333 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, not observable  </td> <td style=\"text-align:left\"> &nbsp;  z[2] contributes to x1 with 100 %<br> </td> </tr>
<tr>
 <td style=\"text-align:center\">       3 </td> <td style=\"text-align:left\"> &nbsp;    2.9891e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3346 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[3] contributes to x2 with 51.9 %<br>&nbsp;  z[3] contributes to x1 with 23.9 % </td> </tr>
<tr>
 <td style=\"text-align:center\">       4 </td> <td style=\"text-align:left\"> &nbsp;    5.5825e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0179 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[4] contributes to x3 with 48.4 %<br>&nbsp;  z[4] contributes to x5 with 32.5 % </td> </tr>
</table>

<p>
<b>Conjugated complex pairs of eigenvalues</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin:20px 0 20px 20px;\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> eigenvalue </td><td> freq. [Hz] </td> <td> damping </td><td> characteristics </td>  <td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:left\">     5/6 </td> <td style=\"text-align:left\"> &nbsp;    1.0299e+000 &plusmn;  6.5528e+000j </td> <td style=\"text-align:left\"> &nbsp;    1.0557 </td> <td style=\"text-align:left\"> &nbsp;   -0.1553 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[    5/6] contribute to x6 with 35.9 %<br>&nbsp;  z[    5/6] contribute to x2 with 20.6 % </td> </tr>
</table>

<p>
In the table above, the column <b>contribution to states</b> lists for each eigenvalue the states
to which thecorresponding modal state contributes most. This information is based on the
two largest absolute values of the corresponding right eigenvector (if the second large value
is less than 5&nbsp;% of the largest contribution, it is not shown).
</p>

<p>
In the next table, for each state in the column <b>correlation to modal states</b>, the modal
states which contribute most to the coresponding state are summarized, i.e. the state is mostly composed of these modal states
This information is based on the two largest absolute values of row i of the
eigenvector matrix that is associated with eigenvalue i (if the second large value
is less than 5&nbsp;% of the largest contribution, it is not shown). This only holds
if the modal states are in the same order of magnitude. Otherwise, the modal states
listed in the last column might be not the most relevant one.
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin:20px 0 20px 20px;\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> state </td> <td> composition </td> <td> eigenvalue #</td> <td> freq. [Hz] </td> <td> damping </td> <td> T [s] </td></tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x1 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  42.5% by z[2] <br> &nbsp;  is composed of  35.4% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 2<br> &nbsp; 5/6 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    1.0557 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp;    0.0201<br> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x2 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  44.2% by z[3] <br> &nbsp;  is composed of  43.7% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 3<br> &nbsp; 5/6 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    1.0557 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp;    0.3333<br> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x3 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  36.9% by z[1] <br> &nbsp;  is composed of  36.3% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 1<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;    0.3346<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x4 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  88.9% by z[5/6] <br> &nbsp;  is composed of   9.8% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 5/6<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp;    1.0557<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;   -0.1553<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x5 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  45.3% by z[1] <br> &nbsp;  is composed of  44.1% by z[4] </td> <td style=\"text-align:center\"> &nbsp; 1<br> &nbsp; 4 </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; ---<br> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp;    0.0000<br> &nbsp;    0.0179 </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp; x6 </td> <td style=\"text-align:left\"> &nbsp;  is composed of  95.7% by z[5/6] </td> <td style=\"text-align:center\"> &nbsp; 5/6          </td> <td style=\"text-align:center\"> &nbsp;    1.0557          </td> <td style=\"text-align:center\"> &nbsp;   -0.1553 </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
</table>

<p>
<b>Invariant zeros</b>
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; margin:20px 0 20px 20px;\" cellpadding=\"3\" border=\"1\" cellspacing=\"0\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> invariant zero </td><td> Time constant [s] </td> <td> freq. [Hz] </td> <td> damping </td></tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;       1 </td> <td> &nbsp;   -5.4983e+001 </td> <td> &nbsp;    0.0182 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;       2 </td> <td> &nbsp;   -3.0000e+000 </td> <td> &nbsp;    0.3333 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;     3/4 </td> <td style=\"text-align:left\"> &nbsp;    3.2417e+000 &plusmn;  5.6548e+000j </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:left\"> &nbsp;    1.0374 </td> <td style=\"text-align:left\"> &nbsp;   -0.4973 </td> </tr>
</table>
</html>", revisions="<html>
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
    end analysis2;
    annotation (Documentation(info="<html>
<p>
This package collects functions used for common analyses on a state space system
represented by a StateSpace record.
</p>
</html>"));
  end Analysis;

  encapsulated package Design
    "Package of functions to design state space controllers and observers"
    import Modelica;
    extends Modelica.Icons.Package;
    encapsulated function assignPolesSI
      "Pole placement for single input systems using Ackermann's formula."

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Complex p[size(ss.A, 1)]
        "Vector of desired poles";
      output Real k[size(ss.A, 1)] "Feedback gain matrix";

    protected
      Real cm[size(ss.B, 1), size(ss.A, 1)*size(ss.B, 2)];
      Modelica_LinearSystems2.Math.Polynomial poly;
      Real Y[size(ss.A, 1), size(ss.A, 2)];
      Real X[:, :];
      Complex p_actual[size(p, 1)];
      Complex p_sorted[size(p, 1)];
      Real poleError;
      Complex smaller;
    algorithm
      assert(size(ss.B, 2) == 1, "System must be SI but has " + String(size(ss.B,
        2)) + " inputs");
      cm := StateSpace.Analysis.controllabilityMatrix(ss);
      assert(Modelica.Math.Matrices.rank(cm) == size(cm, 1) or
        Modelica.Math.Matrices.rank(cm) == size(cm, 2),
        "Controllability matrix has not full rank. System is not controllable!");
      poly := Modelica_LinearSystems2.Math.Polynomial(p);
      Y := Modelica_LinearSystems2.Math.Polynomial.evaluateMatrix(poly, ss.A);
      X := zeros(size(cm, 2), size(Y, 2));
      for j in 1:size(Y, 2) loop
        X[:, j] := Modelica.Math.Matrices.leastSquares(cm, vector(Y[:, j]));
      end for;
      k := X[size(ss.A, 1), :];

      // Check results
      // sort p
      p_sorted := p;
      for i1 in 1:size(p_sorted, 1) loop
        for i2 in (1 + i1):size(p_sorted, 1) loop
          if Modelica.ComplexMath.'abs'(p_sorted[i1]) >
            Modelica.ComplexMath.'abs'(p_sorted[i2]) then
            smaller := p_sorted[i2];
            p_sorted[i2] := p_sorted[i1];
            p_sorted[i1] := smaller;
          end if;
        end for;
      end for;

      p_actual :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(
        ss.A - ss.B*transpose(matrix(k)));
      // sort actual eigenvalues
      for i1 in 1:size(p_actual, 1) loop
        for i2 in (1 + i1):size(p_actual, 1) loop
          if Modelica.ComplexMath.'abs'(p_actual[i1]) >
            Modelica.ComplexMath.'abs'(p_actual[i2]) then
            smaller := p_actual[i2];
            p_actual[i2] := p_actual[i1];
            p_actual[i1] := smaller;
          end if;
        end for;
      end for;

      // check for poles that have an error of more than 10%
      for i in 1:size(p_sorted, 1) loop
        if (Modelica.ComplexMath.'abs'(p_sorted[i]) <> 0) then
          poleError :=Modelica.ComplexMath.'abs'(p_sorted[i] - p_actual[i])/
            Modelica.ComplexMath.'abs'(p_sorted[i]);

          if poleError > 0.1 then
            Modelica.Utilities.Streams.print("Warning: Pole location of pole "
               + String(p_sorted[i]) + " has an error of " + String(100*
              poleError) + "%. (Is " + String(p_actual[i]) + ")");

          end if;
        end if;
      end for;

      annotation (Documentation(revisions="<html>
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
    end assignPolesSI;

    encapsulated function assignPolesMI
      "Pole assignment design algorithm for multi input systems"

      import Modelica;
      // import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;

      input StateSpace ss "State space system";

      input Complex gamma[:]=fill(Complex(0), 0) "Designed Poles";
      //  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
      input Real alpha=-1e10
        "Maximum admissible value for real parts of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
      input Real tolerance=Modelica.Math.Matrices.norm(ss.A, 1)*1e-12
        "Tolerance to be used in determining the controllability of (A,B)";
      input Boolean calculateEigenvectors=false
        "Calculate the eigenvectors X of the closed loop system when true";

      output Real K[size(ss.B, 2), size(ss.A, 1)]
        "State feedback matrix assigning the desired poles";
      output Real S[:, :] "Closed loop System matrix";
      output Complex po[size(ss.A, 1)] "poles of the closed loop system";
      output Integer nfp
        "number of eigenvalues that are not modified with respect to alpha";
      output Integer nap "number of assigned eigenvalues";
      output Integer nup "number of uncontrollable eigenvalues";
      output Complex X[size(ss.A, 1), size(ss.A, 1)]
        "eigenvectors of the closed loop system";

    protected
      Real A_rsf[size(ss.A, 1), size(ss.A, 2)];
      Real B_rsf[size(ss.B, 1), size(ss.B, 2)];
      Real Q[size(ss.A, 1), size(ss.A, 1)];
      Real Ks1[:, :];
      Real Ks2[:, :];
      Real Q2[:, :];
      Real A_rsf_1[:, :];
      Real Q1[:, :];
      Boolean select[:];
      Boolean rselectA[:];
      Real Z[:, :] "orthogonal transformation matrix";
      Real ZT[:, :] "orthogonal transformation matrix";
      Complex pf[:];
      Complex gammaReordered[:]=gamma;
      Integer info;
      Real wr[size(gamma, 1)];
      Real wi[size(gamma, 1)];
      Boolean imag=false;
      Integer i;
      Integer ii;
      Integer iii;
      Integer counter;
      Integer counter2;
      Integer n=size(ss.A, 1);
      Integer nccA "number of conjugated complex pole pairs of openloop system";
      Integer nccg "number of conjugated complex pole pairs of gamma";
      Integer rpg "number of real poles in gamma";
      Integer rpA "number of real poles of open loop system";
      Integer ncc "Min(nccA, nccg)";
      Integer rp "Min(rpg, rpA)";
      Integer ng=size(gamma, 1);
      Integer nr "Differenz between rpA and rpg; Sign(rpA-rpg)*(rpA-rpg)";

      Real alphaReal[size(ss.A, 1)]
        "Real part of eigenvalue=alphaReal+i*alphaImag";
      Real alphaImag[size(ss.A, 1)]
        "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

      Complex SS[:, :];
      Complex Xj[:, :];
      Complex h;

      Real dist;
      Real evImag;

    algorithm
      assert(size(gamma, 1) <= size(ss.A, 1),
        "At most n (order of ss) eigenvalues can be assigned");

      /* Extraction of Poles (Variable conversation) and pole sequence check */
      for i in 1:size(gamma, 1) loop
        wr[i] := gamma[i].re;
        wi[i] := gamma[i].im;
        if imag then
          assert(wi[i - 1] == -wi[i] and wr[i - 1] == wr[i],
            "Poles are in wrong sequence");
          imag := false;
        elseif wi[i] <> 0 then
          imag := true;
        end if;
      end for;

      // put matrix ss.A to real Schur form A <- QAQ' and compute B <- QB
      (A_rsf,Z,alphaReal,alphaImag) := Matrices.rsf2(ss.A);
      ZT := transpose(Z);

      // reorder real Schur form according to alpha
      (A_rsf,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSFc(
            A_rsf,
            identity(size(A_rsf, 1)),
            alphaReal,
            alphaImag,
            alpha);
      ZT := transpose(Z)*ZT;
      B_rsf := ZT*ss.B;

      // determine number of poles not to be assigned according to alpha
      nfp := 0;
      for i in 1:n loop
        if alphaReal[i] < alpha then
          nfp := nfp + 1;
        end if;
      end for;
      nap := n - nfp;

      assert(size(gamma, 1) >= nap, String(nap) +
        " poles should be modified, therefore gamma should contain at at least "
         + String(nap) + " assigned eigenvalues");

      // second reorder (reorderRSF3) according to conjugated complex pairs in A and p
      // count numbre of conjugated complex pole pairs = max(number_ccpp(eig(A), number_ccpp(gamma))
      nccA := 0;
      //mark the real poles of original system
      rselectA := fill(true, nap);
      ii := 1;
      for i in nfp + 1:n loop
        if abs(alphaImag[i]) > 0 then
          nccA := nccA + 1;
        else
          rselectA[ii] := false;

        end if;
        ii := ii + 1;
      end for;
      rpA := n - nccA;
      nccA := div(nccA, 2);

      // reorder gamma and A_rsf
      (gammaReordered,rpg) := Modelica_LinearSystems2.Internal.reorderZeros(
        gamma);
      gammaReordered :=Modelica.ComplexMath.Vectors.reverse(gammaReordered);
      nccg := div(size(gammaReordered, 1) - rpg, 2);
      ncc := min(nccA, nccg);
      rp := min(rpA, rpg);
      if nccA > 0 then
        (A_rsf[nfp + 1:n, nfp + 1:n],Q2) := Matrices.LAPACK.dtrsen(
              "E",
              "V",
              rselectA,
              A_rsf[nfp + 1:n, nfp + 1:n],
              identity(n - nfp));
        //The Schur vector matrix is identity, since A_rsf already has Schur form

        A_rsf[1:nfp, nfp + 1:n] := A_rsf[1:nfp, nfp + 1:n]*Q2;
        B_rsf[nfp + 1:n, :] := transpose(Q2)*B_rsf[nfp + 1:n, :];
        ZT[nfp + 1:n, :] := transpose(Q2)*ZT[nfp + 1:n, :];
      end if;

      // main algorithm
      K := zeros(size(ss.B, 2), size(ss.A, 1));
      counter := nfp + 1;
      counter2 := 1;

      for i in 1:rp loop
        // 1x1 blocks; real system pole and real assigned poles; take the next eigenvalue in the
        // diagonal of the Schur form and search the nearest pole in the set of the real poles to assign
        dist := Modelica.Constants.inf;
        for ii in i:rpg loop
          // looking for nearest pole and reorder gamma
          if abs(A_rsf[n, n] - gammaReordered[ng - ii + 1].re) < dist then
            iii := ng - ii + 1;
            dist := abs(A_rsf[n, n] - gammaReordered[ng - ii + 1].re);
          end if;
        end for;
        h := gammaReordered[ng - i + 1];
        gammaReordered[ng - i + 1] := gammaReordered[iii];
        gammaReordered[iii] := h;

        Ks1 := StateSpace.Internal.assignOneOrTwoPoles(
              matrix(A_rsf[n, n]),
              transpose(matrix(B_rsf[n, :])),
              {gammaReordered[ng - i + 1]},
              tolerance);
        K := K + [zeros(size(Ks1, 1), size(K, 2) - 1), Ks1]*ZT;
        A_rsf := A_rsf - B_rsf*[zeros(size(Ks1, 1), size(K, 2) - 1), Ks1];
        select := fill(false, n - counter + 1);
        select[n - counter + 1] := true;

        (A_rsf[counter:n, counter:n],Q1) := Matrices.LAPACK.dtrsen(
              "E",
              "V",
              select,
              A_rsf[counter:n, counter:n],
              identity(n - counter + 1));
        //The Schur vector matrix is identity, since A_rsf already has Schur form

        A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q1;
        B_rsf[counter:n, :] := transpose(Q1)*B_rsf[counter:n, :];
        ZT[counter:n, :] := transpose(Q1)*ZT[counter:n, :];
        counter := counter + 1;
        counter2 := counter2 + 1;
      end for;

      if counter2 < rpg and counter2 > rpA then
        //System has less real eigenvalues than real assigned poles
        for i in 1:div(rpg - rpA, 2) loop
          // 2x2 blocks; complex pair of system poles and 2 real assigned poles; take the next complex pair
          // (Schur bump) in the diagonal of the Schur form and search the two nearest poles in the set of the
          // remaining real assigned poles
          dist := Modelica.Constants.inf;
          evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);
          //positive imaginary part of the complex system pole pair
          for ii in 2*(i - 1) + 1:2:rpg - rpA loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag
                 < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) +
                evImag;
            end if;
          end for;
          h := gammaReordered[ng - rp - 2*(i - 1)];
          gammaReordered[ng - rp - 2*(i - 1)] := gammaReordered[iii];
          gammaReordered[iii] := h;
          dist := Modelica.Constants.inf;
          for ii in 2*(i - 1) + 1:2:rpg - rpA loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag
                 < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) +
                evImag;
            end if;
          end for;
          h := gammaReordered[ng - rp - 2*i + 1];
          gammaReordered[ng - rp - 2*i + 1] := gammaReordered[iii];
          gammaReordered[iii] := h;

          Ks2 := StateSpace.Internal.assignOneOrTwoPoles(
                A_rsf[n - 1:n, n - 1:n],
                matrix(B_rsf[n - 1:n, :]),
                gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)],
                tolerance);

          K := K + [zeros(size(Ks2, 1), size(K, 2) - 2), Ks2]*ZT;
          A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2), Ks2];
          select := fill(false, n - counter + 1);
          select[n - counter:n - counter + 1] := {true,true};

          (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
                "E",
                "V",
                select,
                A_rsf[counter:n, counter:n],
                identity(n - counter + 1));
          //The Schur vector matrix is identity, since A_rsf already has Schur form

          A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
          B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
          ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
          counter := counter + 2;
          counter2 := counter2 + 2;
        end for;
      end if;

      if counter2 > rpg and counter2 < rpA then
        //System has more real eigenvalues than real assigned poles
        for i in 1:div(rpA - rpg, 2) loop
          // 2x2 blocks; 2 real system poles and a pair of complex assigned poles; take the next two real
          // eigenvalues in the diagonal of the Schur form and search the complex pole pair of the assigned poles
          // which is nearest to the two real poles
          dist := Modelica.Constants.inf;
          for ii in 2*(i - 1) + 1:2:rpA - rpg loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + abs(
                gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[n - 1, n - 1]
                 - gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng
                 - rp - ii + 1].im) < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) +
                abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[n - 1, n
                 - 1] - gammaReordered[ng - rp - ii + 1].re) + abs(
                gammaReordered[ng - rp - ii + 1].im);
            end if;
          end for;
          h := gammaReordered[ng - rp - 2*(i - 1)];
          gammaReordered[ng - rp - 2*(i - 1)] := gammaReordered[iii];
          gammaReordered[iii] := h;
          h := gammaReordered[ng - rp - 2*i + 1];
          gammaReordered[ng - rp - 2*i + 1] := gammaReordered[iii - 1];
          gammaReordered[iii - 1] := h;

          Ks2 := StateSpace.Internal.assignOneOrTwoPoles(
                A_rsf[n - 1:n, n - 1:n],
                matrix(B_rsf[n - 1:n, :]),
                gammaReordered[ng - rp - 2*i + 1:ng - rp - 2*(i - 1)],
                tolerance);

          K := K + [zeros(size(Ks2, 1), size(K, 2) - 2), Ks2]*ZT;
          A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2), Ks2];
          select := fill(false, n - counter + 1);
          select[n - counter:n - counter + 1] := {true,true};

          (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
                "E",
                "V",
                select,
                A_rsf[counter:n, counter:n],
                identity(n - counter + 1));
          //The Schur vector matrix is identity, since A_rsf already has Schur form

          A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
          B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
          ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
          counter := counter + 2;
          counter2 := counter2 + 2;
          //      Modelica.Utilities.Streams.print("counter2Case3 = " + String(counter2));
        end for;
      end if;

      for i in 1:ncc loop
        // 2x2 blocks; 2 complex system poles and two complex assigned poles; take the next complex
        // system pole pair (next Schur bump) in the diagonal of the Schur form and search the complex
        //  assigned pole pair which is nearest
        dist := Modelica.Constants.inf;
        evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);
        //positive imaginary part of the complex system pole pair
        for ii in 2*(i - 1) + 1:2:2*ncc loop
          if abs(A_rsf[n, n] - gammaReordered[2*ncc - ii + 1].re) + abs(evImag
               - abs(gammaReordered[2*ncc - ii + 1].im)) < dist then
            iii := 2*ncc - ii + 1;
            dist := abs(A_rsf[n, n] - gammaReordered[2*ncc - ii + 1].re) + abs(
              evImag - abs(gammaReordered[2*ncc - ii + 1].im));
          end if;
        end for;
        h := gammaReordered[2*ncc - 2*(i - 1)];
        gammaReordered[2*ncc - 2*(i - 1)] := gammaReordered[iii];
        gammaReordered[iii] := h;
        h := gammaReordered[2*ncc - 2*i + 1];
        gammaReordered[2*ncc - 2*i + 1] := gammaReordered[iii - 1];
        gammaReordered[iii - 1] := h;

        Ks2 := StateSpace.Internal.assignOneOrTwoPoles(
              A_rsf[n - 1:n, n - 1:n],
              matrix(B_rsf[n - 1:n, :]),
              gammaReordered[2*ncc - 2*i + 1:2*ncc - 2*(i - 1)],
              tolerance);
        K := K + [zeros(size(Ks2, 1), size(K, 2) - 2), Ks2]*ZT;
        A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2), Ks2];
        select := fill(false, n - counter + 1);
        select[n - counter:n - counter + 1] := {true,true};

        (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
              "E",
              "V",
              select,
              A_rsf[counter:n, counter:n],
              identity(n - counter + 1));
        //The Schur vector matrix is identity, since A_rsf already has Schur form

        A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
        B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
        ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
        counter := counter + 2;
        counter2 := counter2 + 2;
      end for;

      S := ss.A - ss.B*K;
      po :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(S);

      if calculateEigenvectors then
        //     X := fill(Complex(0), n, n);
        //     for i in 1:n loop
        //       SS := Complex(1)*S;
        //       for ii in 1:n loop
        //         SS[ii, ii] := SS[ii, ii] - po[i];
        //       end for;
        //       Xj := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_nullspace(
        //                                  SS);
        //       for ii in 1:n loop
        //         X[ii, i] := Xj[ii, 1];
        //       end for;
        //     end for;
        //      Modelica_LinearSystems2.Math.Complex.Matrices.print(X,6,"X1");
        X :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenVectors(S);
        //      Modelica_LinearSystems2.Math.Complex.Matrices.print(X,6,"X2");

      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(K, S, po, nfp, nap, nup) = StateSpace.Design.<b>assignPolesMI</b>(ss, gamma, np, tol, calculateEigenvectors)
</pre></blockquote>

<h4>Description</h4>
<p>
The purpose of this function is to determine the state feedback matrix <b>K</b> for a
given time invariant multi input state system (<b>A</b>,<b>B</b>) such that the
closed-loop state matrix <b>A</b>-<b>B</b>*<b>K</b> has specified eigenvalues. The
feedback matrix <b>K</b> is calculated by factorization following [1]. The algorithm
modifies the eigenvalues sequentially and also allows partial eigenvalue assignment.
</p>
<p>
At the beginning of the algorithm, the feedback matrix <b>K</b> is set to zero (<b>K</b> = <b>0</b>) and the matrix <b>A</b> is
reduced to an ordered real Schur form by separating its spectrum in two parts
</p>
<blockquote><pre>
             | <b>F</b>1  <b>F</b>3|
<b>F</b> = <b>Q</b>*<b>A</b>*<b>Q</b>' = |       |
             | <b>0</b>   <b>F</b>2|
</pre></blockquote>
<p>
in such a way, that <b>F</b>1 contains the eigenvalues that will be
retained and <b>F</b>3 contains the eigenvalues going to be modified. On the suggestion
of [1] the eigenvalues <i>evr</i> to be retained are chosen as
</p>
<blockquote><pre>
evr = {s in C: Re(s) &lt; -alpha, alpha &gt;= 0}
</pre> </blockquote>
<p>
but other specification are conceivable of course.
</p>
<p>
Let
</p>
<blockquote><pre>
<b>G</b> = [<b>G</b>1;<b>G</b>2] = <b>Q</b>*<b>B</b>
</pre> </blockquote>
<p>
with an appropriate partition according to <b>F</b>2. (<b>F</b>2, <b>G</b>2) has to be
controllable.
</p>
<p>
If the feedback matrix <b>K</b> is taken in a form
</p>
<blockquote><pre>
<b>K</b> = [0, <b>K</b>2]
</pre></blockquote>
<p>
the special structure of <b>F</b> and <b>K</b> results in a closed loop state
matrix
</p>
<blockquote><pre>
          |<b>F</b>1 <b>F</b>3 - <b>G</b>1*<b>K</b>2|
<b>F</b> - <b>G</b>*<b>K</b> = |             |
          |0  <b>F</b>2 - <b>G</b>2*<b>K</b>2|
</pre></blockquote>
<p>
with only the eigenvalues of <b>F</b>2 are modified. This approach to modify
separated eigenvalues is used to sequentially shift one real eigenvalue ore two
complex conjugated eigenvalues stepwise until all assigned eigenvalues are placed.
Therefore, at each step i always the (two) lower right eigenvalue(s) are modified by an
appropriate feedback matrix <b>K</b>i. The matrix <b>F</b> - <b>G</b>*<b>K</b>i remains in real Schur form. The
assigned eigenvalue(s) is (are) then moved to another diagonal position of the real Schur
form using reordering techniques <b>F</b> &lt; -- <b>Q</b>i*<b>F</b>*<b>Q</b>i'  and a new block is transferred to the
lower right diagonal position. The transformations are accumulated in <b>Q</b>i and are also
applicated to the matrices
</p>
<blockquote><pre>
<b>G</b> &lt; - <b>Q</b>i*<b>G</b> <b>Q</b> &lt; - <b>Q</b>i*<b>Q</b>
</pre></blockquote>
<p>
The eigenvalue(s) to be assigned at  each step is (are) chosen such that the norm of each <b>K</b>i is minimized [1].
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1,1; 0,-2],
    B=[0; 1],
    C=[1,0; 0,1],
    D=[0; 0]);

  Complex p[:]={Complex(-3,0),Complex(-4,0)};

<b>algorithm</b>
  (K, S, newPoles) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, p);

 // K = [6.0, 4.0]
 // S = [-1.0, 1.0; -6.0, -6.0]
 // newPoles = {-3, -4}
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga A. (1981):</dt>
<dd> <b>A Schur method for pole assignment</b>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
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
    end assignPolesMI;

    function kalmanFilter "Design of a Kalman estimator matrix"
      import Complex;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "Time-continuous system in state space form";
      input Real Q[size(ss.A, 1), size(ss.A, 1)]
        "Covariance Matrix of state noise (n x n), n number of states";
      input Real R[size(ss.C, 1), size(ss.C, 1)]
        "Covariance Matrix of output noise (m x m), m number of inputs";

    public
      output Real L[:, :] "Kalman filter matrix";
      output StateSpace kss(
        redeclare Real A[size(ss.A, 1), size(ss.A, 1)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2) + size(ss.C, 1)],
        redeclare Real C[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real D[size(ss.A, 1), size(ss.B, 2) + size(ss.C, 1)])
        "kalman system";

    protected
      Real AR[:, :]=transpose(ss.A);
      Real BR[:, :]=transpose(ss.C);
      Real CR[:, :]=zeros(1, size(AR, 1));
      Real DR[:, :]=zeros(1, size(BR, 2));
      StateSpace rss=StateSpace(
              AR,
              BR,
              CR,
              DR) "System to calculate the Kalman estimator with lqr algorithm";

      Real AK[size(ss.A, 1), size(ss.A, 1)];
      Real BK[size(ss.B, 1), size(ss.B, 2) + size(ss.C, 1)];
      Real CK[size(ss.A, 1), size(ss.A, 2)];
      Real DK[size(ss.A, 1), size(ss.B, 2) + size(ss.C, 1)];
      // matrices of the kalman system kss

    algorithm
      (L) := StateSpace.Design.lqr(
            rss,
            Q,
            R);
      L := transpose(L);

      AK := ss.A - L*ss.C;
      BK[:, 1:size(ss.B, 2)] := ss.B - L*ss.D;
      BK[:, size(ss.B, 2) + 1:size(BK, 2)] := L;
      CK := identity(size(ss.A, 1));
      DK := zeros(size(ss.A, 1), size(ss.B, 2) + size(ss.C, 1));

      kss := StateSpace(
            AK,
            BK,
            CK,
            DK);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(L, kss) = StateSpace.Design.<b>kalmanFilter</b>(ss, Q, R)
</pre></blockquote>

<h4>Description</h4>
<p>
This functions designs the kalman-bucy filter, that reconstructs
plant states and system output without noise.
As input it uses the plant input and output.
</p>
<p>
Noise affects the plant states via q(t)
</p>
<blockquote>d
x/dt = Ax + Bu + q(t)
</blockquote>
<p>
The plant output is affected by r(t)
</p>
<blockquote>
y = Cx + Du + r(t)
</blockquote>
<p>
The covariance matrices of q and r have to be given via Q and R, respectively.
</p>
<p>
The filter uses an observer that tries to reconstruct the original behaviour. Its states and outputs are trailed with a hat (^)<br>
The observer is controlled by feedback of the output difference y - y^ (y^= Cx^+ Du)
over a Matrix L, such that x^
</p>
<blockquote>
dx^/dt = (A - LC) x^ + (B - LD)u + Ly
</blockquote>
<p>
follows the plant state x. L is designed to minimize noise in states and inputs.
L is calculated from a Riccati Equation
</p>
<blockquote><pre>
    -1
SC'R  CS - SA' - AS - Q = 0
        -1
L = SC'R
</pre></blockquote>

<p>
The representation of the estimation model would be as follows:
</p>
<blockquote><pre>
.                            |u|
x^  = [A-LC] x^ + [B-LD , L] | |
                             |y|
                    |u|
y^ = [C] x^ + [D 0] | |
                    |y|
</pre></blockquote>
<p>
Since the controller approach was made to provide the estimated states, the representation of the ooutput kss is such that
</p>
<blockquote><pre>
y^ = x^
</pre></blockquote>
<p>
i.e., kss:
</p>
<blockquote><pre>
.                            |u|
x^  = [A-LC] x^ + [B-LD , L] | |
                             |y|
            |u|
y^ = Ix^ + 0| |
            |y|
</pre></blockquote>
<p>
i.e.
</p>
<blockquote><pre>
C^ = I,   D^ = 0   with appropriate sizes, i.e. size(C^) = {nx,nx},  size(D^) = {nx, nu+ny}.
</pre></blockquote>

<p>
Since the calculation of a Kalman filter is the dual problem of lqr calculation function
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Design.lqr\">Modelica_LinearSystems2.StateSpace.Design.lqr</a>
is used to solve the Riccati euation.<br>
The algebraic Riccati equation is solved by using the Schur algorithm
<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.TransferFunction;
  StateSpace css = StateSpace(TransferFunction.constructor({1},{1,2,3,4}));
  Real Q[:,:] = identity(3);
  Real R[:,:] = [1];

  Real L[:,:];
  StateSpace kss(
            redeclare Real A[size(css.A,1),size(css.A,1)],
            redeclare Real B[size(css.B,1),size(css.B,2)+size(css.C,1)],
            redeclare Real C[size(css.C,1)+size(css.A,1),size(css.C,2)],
            redeclare Real D[size(css.C,1)+size(css.A,1),size(css.B,2)+size(css.C,1)]);

<b>algorithm</b>
  ( L, kss) := StateSpace.Design.kalmanFilter(css, Q, R);
//  L = [0.9928;
       (-0.0072);
       (-1.4868)]
// kss = StateSpace(
          A = [-0.993,     1,     0;
                0.007,     0,     1;
               -2.513,     3,    -2],

          B = [ 0,         0.993;
                0,        -0.0072;
                1,        -1.4868],

          C = [ 1,     0,     0;
                0,     1,     0;
                0,     0,     1],

          D = [ 0,     0;
                0,     0;
                0,     0])
</pre></blockquote>
</html>", revisions="<html>
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
    end kalmanFilter;

    encapsulated function lqr "LQR design algorithm"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math;

      input StateSpace ss "Open loop system in state space form";
      input Real Q[size(ss.A, 1), size(ss.A, 2)]=identity(size(ss.A, 1))
        "State weighting matrix";
      input Real R[size(ss.B, 2), size(ss.B, 2)]=identity(size(ss.B, 2))
        "Input weighting matrix";
    protected
      Boolean iscontinuousSystem=true;
    public
      output Real K[size(ss.B, 2), size(ss.A, 1)] "Feedback gain matrix";
      output StateSpace sslqr(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]) "Closed loop system";

      output Real S[size(ss.A, 1), size(ss.A, 1)] "Solution of the Riccati equation";
      output Complex ev[:] "Eigenvalues of the system";

    algorithm
      if min(size(ss.A, 1), size(ss.A, 2)) > 0 then
        assert(StateSpace.Analysis.isControllable(ss),
          "System in function \"Modelica_LinearSystems2.StateSpace.Design.lqr\" has to be controllable");
        if iscontinuousSystem then
          (S,ev) := Math.Matrices.care(
                ss.A,
                ss.B,
                R,
                Q);
          K := Math.Matrices.solve2(R, transpose(ss.B)*S);
        else
          (S,ev) := Math.Matrices.dare(
                ss.A,
                ss.B,
                R,
                Q);
          K := Math.Matrices.solve2(R + transpose(ss.B)*S*ss.B, transpose(ss.B)
            *S*ss.A);
        end if;

        sslqr.A := ss.A - ss.B*K;
        sslqr.B := ss.B;
        sslqr.C := ss.C - ss.D*K;
        sslqr.D := ss.D;

      else
        K := fill(
              0,
              size(ss.B, 2),
              size(ss.A, 1));
        ev := fill(Complex(0), 0);
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(K, sslqr, X, ev) = StateSpace.<b>lqr</b>(ss, Q, R, true)
</pre></blockquote>

<h4>Description</h4>
<p>
The optimal and stabilizing gain matrix <b>K</b> for a state-feedback law <b>u</b> = -<b>K</b>*<b>x</b>
is designed such that the cost function
</p>
<blockquote><pre>
J = Integral {<b>x</b>'*<b>Q</b>*<b>x</b> + <b>u</b>'*<b>R</b>*<b>u</b>} dt
</pre></blockquote>
<p>
of the continuous time case or
</p>
<blockquote><pre>
Jd = Sum {<b>x</b>'k*<b>Q</b>*<b>x</b>k + <b>u</b>'k*<b>R</b>*<b>u</b>k}
</pre></blockquote>
<p>
of the discrete time case is minimized. The cases are chosen by the input <b>iscontinuousSystem</b> This is done by solving
the continuous-time algebraic Riccati equation (CARE)
</p>
<blockquote><pre>
<b>Q</b> + <b>A</b>'*<b>X</b> + <b>X</b>*<b>A</b> - <b>X</b>*<b>B</b>*<b>R</b><sup>-1</sup>*<b>B</b>'*<b>X</b> = <b>0</b>
</pre></blockquote>
<p>
or the discrete-time algebraic Riccati equation (DARE)
</p>
<blockquote><pre>
<b>X</b> - <b>A</b>'*<b>X</b>*<b>A</b> + <b>A</b>'*<b>X</b>*<b>B</b>*(<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)<sup>-1</sup>*<b>B</b>'*<b>X</b>*<b>A</b> - <b>Q</b> = <b>0</b>
</pre>
</blockquote>
<p>
for <b>X</b> using the Schur vector approach. See <a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a> and <a href=\"Modelica://Modelica_LinearSystems2.Math.Matrices.dare\">dare</a> respectively for more details.
</p>
<p>
The gain matrix <b>K</b> of the continuous-time case is calculated from
</p>
<blockquote><pre>
<b>K</b> = <b>R</b><sup>-1</sup>*<b>B</b>'*<b>X</b>
</pre></blockquote>
<p>
or from
</p>
<blockquote><pre>
<b>K</b> = (<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)<sup>-1</sup>*<b>B</b>'*<b>X</b>*<b>A</b>
</pre></blockquote>
<p>
for the discrete-time case.
The output state space system sslqr represents the closed loop system
</p>
<blockquote><pre>
  .
  <b>x</b> = [<b>A</b> - <b>BK</b>] <b>x</b> + <b>Bu</b>

  <b>y</b> = [<b>C</b> - <b>DK</b>] <b>x</b> + <b>Du</b>

</pre></blockquote>
<p>
The output S is the solution of the Riccati equation
</p>
<p>
The eigenvalues of the closed loop system <b>A</b> - <b>B</b>*<b>K</b> are computed as complex output ev.
</p>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss=StateSpace(
    A=[0, 1, 0, 0; 0, 0, 39.2, 0; 0, 0, 0, 1; 0, 0, 49, 0],
    B=[0; 1; 0; 1],
    C=[1, 0, 0, 0],
    D=[0]);
  Real Q[:,:]=identity(4);
  Real R[:,:]=identity(1);
  Real K[size(ss.B, 2),size(ss.A, 1)];

<b>algorithm</b>
  K := StateSpace.Design.lqr(ss, Q, R);

// K = [-1, -3.63271, 108.763, 18.3815]
</pre></blockquote>
</html>", revisions="<html>
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
    end lqr;

    encapsulated function lqg "LQG design algorithm"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math;

      input StateSpace ss "Open loop system in state space form";
      input Real Q[size(ss.A, 1), size(ss.A, 2)]=identity(size(ss.A, 1))
        "State weighting matrix";
      input Real R[size(ss.B, 2), size(ss.B, 2)]=identity(size(ss.B, 2))
        "Input weighting matrix";
      input Real V[size(ss.C, 1), size(ss.C, 1)]=identity(size(ss.C, 1))
        "Covariance output noise matrix";
      input Real W[size(ss.A, 1), size(ss.A, 1)]=identity(size(ss.A, 1))
        "Covariance state noise matrix";

      input Boolean iscontinuousSystem=true
        "True, if state space system is continuous";

      output Real Kc[size(ss.B, 2), size(ss.A, 1)]
        "Controller feedback gain matrix";
      output Real Kf[size(ss.A, 1), size(ss.C, 1)]
        "Kalman feedback gain matrix";

      output StateSpace sslqg(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.C, 1)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.C, 1), size(ss.C, 1)]) "Closed loop system";

    protected
      Real AR[:, :]=transpose(ss.A);
      Real BR[:, :]=transpose(ss.C);
      Real CR[:, :]=zeros(1, size(AR, 1));
      Real DR[:, :]=zeros(1, size(BR, 2));
      // System for lqr
      StateSpace rss=StateSpace(
              AR,
              BR,
              CR,
              DR);
      Real Sc[size(ss.A, 1), size(ss.A, 1)]
        "Solution of the Riccati equation, controller";
      Real Sf[size(ss.A, 1), size(ss.A, 1)]
        "Solution of the Riccati equation, filter";

    algorithm
      if min(size(ss.A, 1), size(ss.A, 2)) > 0 then
        assert(StateSpace.Analysis.isControllable(ss),
          "System in function \"Modelica_LinearSystems2.StateSpace.Design.lqg\" has to be controllable");

        if iscontinuousSystem then
          (Sc,) := Math.Matrices.care(
                ss.A,
                ss.B,
                R,
                Q);
          Kc := Math.Matrices.solve2(R, transpose(ss.B)*Sc);
        else
          (Sc,) := Math.Matrices.dare(
                ss.A,
                ss.B,
                R,
                Q);
          Kc := Math.Matrices.solve2(R + transpose(ss.B)*Sc*ss.B, transpose(ss.B)
            *Sc*ss.A);
        end if;

        assert(StateSpace.Analysis.isObservable(ss),
          "System in function \"Modelica_LinearSystems2.StateSpace.Design.lqg\" has to be observable");
        if iscontinuousSystem then
          (Sf,) := Math.Matrices.care(
                rss.A,
                rss.B,
                V,
                W);
          Kf := transpose(Math.Matrices.solve2(V, ss.C*Sf));
        else
          (Sf,) := Math.Matrices.dare(
                rss.A,
                rss.B,
                V,
                W);
          Kf := transpose(Math.Matrices.solve2(V + rss.C*Sf*rss.B, rss.C*Sf*rss.A));
        end if;

      else
        Kc := fill(
              0,
              size(ss.B, 2),
              size(ss.A, 1));
        Kf := fill(
              0,
              size(ss.C, 1),
              size(ss.A, 1));

      end if;

      sslqg.A := ss.A - Kf*ss.C - ss.B*Kc + Kf*ss.D*Kc;
      sslqg.B := Kf;
      sslqg.C := ss.C - ss.D*Kc;
      sslqg.D := zeros(size(ss.C, 1), size(ss.C, 1));

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(Kc, Kf, sslqg) = StateSpace.<b>lqg</b>(ss, Q, R, V, W)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates matrices <b>K</b>c and <b>K</b>f for linear quadratic gaussian problem (LQG), i.e. the minimization of the expected value of a cost function in consideration of stochastically disturbed states and outputs of the system
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b><b>x</b> + <b>B</b><b>u</b> + <b>w</b>
     <b>y</b> = <b>C</b><b>x</b> + <b>D</b><b>u</b> + <b>v</b>
</pre></blockquote>
<p>
The noise <b>w</b>(t) and <b>v</b>(t) are supposed to be both white, Gaussian zero-mean, stationary stochastic processes with positive semidefinte covariance matrix <b>W</b>
</p>
<blockquote><pre>
E[<b>w</b>(t)*<b>w</b>'(tau)] = <b>W</b>*delta(t-tau)
</pre></blockquote>
<p>
and positive covariance matrix <b>V</b>
</p>
<blockquote><pre>
E[<b>v</b>(t)*<b>v</b>'(tau)] = <b>V</b>*delta(t-tau).
</pre></blockquote>
<p>
E[s] denotes the expected value of a signal s.
</p>
<p>
The LQG approach combines the deterministic <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Design.lqr\">LQR</a> approach and <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Design.kalmanFilter\">Kalman filter</a> principle to estimate stochastically disturbed systems, such that input <b>u</b>(t) is given by
</p>
<blockquote><pre>
<b>u</b>(t) = -<b>K</b>c<b>x</b>^(t)
</pre></blockquote>
<p>
where <b>K</b>c is a lqr feedback matrix and x^(t) the reconstructed state vector estimated by a Kalman filter.
</p>
<p>
Since, the considered problem is stochastic, the objective function to minimize is an expected value
</p>
<blockquote><pre>
           1      T
J = lim   ---- E[Integral (<b>x</b>'*<b>Q</b>*<b>x</b> + <b>u</b>'*<b>R</b>*<b>u</b>)dt],
  (T->inf) 2T    -T
</pre></blockquote>
<p>
where the weighting matrices <b>Q</b> and <b>R</b> are, respectively, symmetric positive semidefinite and positive definite.
</p>
<p>
The feedback matrix Kc is calculated by
</p>
<blockquote><pre>
<b>K</b>c = <b>R</b><sup>-1</sup>*<b>B</b>'*<b>X</b>c,
</pre></blockquote>
<p>
where <b>X</b>c satisfying the continuous-time algebraic Riccati equation (<a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a>)
</p>
<blockquote><pre>
<b>Q</b> + <b>A</b>'*<b>X</b>c + <b>X</b>c*<b>A</b> - <b>X</b>c*<b>B</b>*<b>R</b><sup>-1</sup>*<b>B</b>'*<b>X</b>c = <b>0</b>.
</pre></blockquote>
<p>
The matrix <b>K</b>f of the filter problem to generate the estimated state vector <b>x</b>^(t) is given by
</p>
<blockquote><pre>
<b>K</b>f = <b>X</b>f*<b>C</b>T*<b>V</b>-1,
</pre></blockquote>
<p>
where <b>X</b>f is satisfying the continuous-time algebraic Riccati equation
</p>
<blockquote><pre>
<b>W</b> + <b>A</b>*<b>X</b>f + <b>X</b>f*<b>A</b>' - <b>X</b>f*<b>C</b>'*<b>V</b><sup>-1</sup>*<b>C</b>*<b>X</b>f = <b>0</b>.
</pre></blockquote>
<p>
The vector <b>x</b>^(t) satisfies the differential equation
</p>
<blockquote><pre>
.
<b>x</b>^(t) = (<b>A</b> - <b>K</b>f<b>C</b>)<b>x</b>^(t) + (<b>B</b> - <b>K</b>f<b>D</b>)<b>u</b>(t) + <b>K</b>f<b>y</b>(t)
</pre></blockquote>
<p>
Combining the equation state feedback and state estimation, the state vector <b>x</b>(t) and the estimated state vector <b>x</b>^(t) are given by
</p>
<blockquote><pre>
 .
|<b>x</b> |   | <b>A</b>         -<b>B</b><b>K</b>c      |  |<b>x</b> |   | <b>I</b>   <b>0</b> |  | <b>w</b> |
|  | = |                     |  |  | + |       |  |   |
|<b>x</b>^|   | <b>K</b>f<b>C</b>   <b>A</b> - <b>B</b><b>K</b>c - <b>K</b>f<b>C</b> |  |<b>x</b>^|   | <b>0</b>  <b>K</b>f |  | <b>v</b> |.
</pre></blockquote>
<p>
Finally, the output sslqg represents the estimated system with <b>y</b>(t), the output of the real system, as the input
</p>
<blockquote><pre>
.
<b>x</b>^ = [<b>A</b> - <b>K</b>f<b>C</b> - <b>B</b><b>K</b>c + <b>K</b>f<b>D</b><b>K</b>c]*<b>x</b>^ + <b>K</b>f*<b>y</b>

<b>y</b>^ = [<b>C</b> - <b>D</b><b>K</b>c] <b>x</b>^
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  StateSpace ss=StateSpace(
    A=[-0.02, 0.005, 2.4,  -32; -0.14,  0.44,  -1.3,  -30; 0,  0.018,  -1.6,  1.2; 0, 0, 1, 0],
    B=[0.14,  -0.12; 0.36, -8.6; 0.35, 0.009; 0, 0],
    C=[0, 1, 0, 0; 0, 0, 0, 57.3],
    D=[0,0; 0,0]);

   Real Q[:,:] = transpose(ss.C)*ss.C \" state weighting matrix\";
   Real R[:,:] = identity(2) \" input weighting matrix\";
   Real V[:,:] = identity(2) \" covariance output noise matrix\";
   Real W[:,:] = ss.B*transpose(ss.B) \" covariance state noise matrix\";
   Real Kc[size(ss.B, 2),size(ss.A, 1)] \"Controller feedback gain matrix\";
   Real Kf[size(ss.A, 1),size(ss.C, 1)] \"Kalman feedback gain matrix\";

<b>algorithm</b>
  (Kc, Kf) := StateSpace.Design.lqg(ss, Q, R, V, W);

// Kc = [-0.0033,     0.04719,      14.6421,        60.8894;
          0.0171,    -1.05154,       0.29273,        3.2468]

// Kf = [0.015,     -0.2405;
         9.066,     -0.1761;
         0.009,      0.2289;
        -0.003,      0.08934]
</pre></blockquote>
</html>", revisions="<html>
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
    end lqg;

  end Design;

  encapsulated package Plot
    "Package of functions to plot state space system responses"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function polesAndZeros
      "Plot poles (i.e. eigenvalues) and/or invariant zeros of a state space system"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "Linear system in state space form"
        annotation (Dialog);
      input Boolean poles=true
        "= true, to plot the poles (i.e. the eigenvalues) of ss"
        annotation (choices(checkBox=true));
      input Boolean zeros=true "= true, to plot the (invariant) zeros of ss "
        annotation (choices(checkBox=true));

      input Boolean print=true
        "= true, to print the selection to the output window"
        annotation (choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            if poles and zeros then
              Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros()
            else if poles then
              Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
                heading="Eigenvalues (x)")
            else
              Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
                heading="Invariant zeros (o)"));
    protected
      Integer nx=size(ss.A, 1);
      Real EigReal[:, 2];
      Real InvZerosReal[:,2];
      Complex invZeros[:];
      Complex eig[:];
      Plot.Records.Curve curves[2];
      Integer i;
      Plot.Records.Diagram diagram2;
      Real Ab[size(ss.A,1),size(ss.A,1)];
      Real Bb[size(ss.A,1),size(ss.B,2)];
      Real Cb[size(ss.C,1),size(ss.A,1)];
    algorithm
      if poles and size(ss.A, 1) > 0 then
        EigReal := Modelica_LinearSystems2.Math.Matrices.eigenValuesAsRealMatrix(ss.A);
      else
        EigReal := fill(0.0,0,2);
      end if;

      if zeros and size(ss.A,1) > 0 and size(ss.B,2) > 0 and size(ss.C,1) > 0 then
        (,Ab,Bb,Cb) :=Modelica_LinearSystems2.Internal.balanceABC(
          ss.A,
          ss.B,
          ss.C);
        InvZerosReal := Modelica_LinearSystems2.StateSpace.Internal.invariantZerosWithRealMatrix(Ab,Bb,Cb,ss.D);
      else
        InvZerosReal :=fill(0.0,  0, 2);
      end if;

      i := 0;
      if poles then
        i := i + 1;
        curves[i] := Plot.Records.Curve(
              x=EigReal[:, 1],
              y=EigReal[:, 2],
              legend="Eigenvalues",
              lineColor={0,0,255},
              autoLine=false,
              linePattern=Plot.Types.LinePattern.None,
              lineSymbol=Plot.Types.PointSymbol.Cross);
      end if;
      if zeros then
        i := i + 1;
        curves[i] := Plot.Records.Curve(
              x=InvZerosReal[:,1],
              y=InvZerosReal[:,2],
              legend="invariant zeros",
              lineColor={255,0,0},
              autoLine=false,
              linePattern=Plot.Types.LinePattern.None,
              lineSymbol=Plot.Types.PointSymbol.Circle);
      end if;

      diagram2 := defaultDiagram;
      diagram2.curve := curves[1:i];
      Plot.diagram(diagram2, device);

      if print then
        if poles then
          eig :=fill(Complex(0), size(EigReal, 1));
          for i in 1:size(eig,1) loop
            eig[i].re :=EigReal[i, 1];
            eig[i].im :=EigReal[i, 2];
          end for;
          Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.printHTML(
            eig, heading="Eigenvalues", name="eigenvalue");
        end if;

        if zeros then
          invZeros :=fill(Complex(0), size(InvZerosReal, 1));
          for i in 1:size(invZeros,1) loop
            invZeros[i].re :=InvZerosReal[i, 1];
            invZeros[i].im :=InvZerosReal[i, 2];
          end for;
          Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.printHTML(
            invZeros, heading="Invariant zeros", name="invariant zero");
        end if;
      end if;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>polesAndZeros</b>(ss);
   or
StateSpace.Plot.<b>polesAndZeros</b>(
  ss,
  poles=true,
  zeros=true,
  print=true,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>());
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a pole-zero-map of the poles and transmission zeros of a state space system.
The poles are the eigenvalues of the system matrix (eigenvalues(ss.A)). The Boolean inputs
\"poles\" and \"zeros\" define what to plot. If Boolean input \"plot = true\", the pole-zero-map
is plotted. If false, only the diagram is generated and returned as output argument.
When the Boolean input print=true, then the results of the analysis are additionally printed
in the textual format.
The records \"defaultDiagram\" and \"device\" allow to set various layout options and the
size and location of the diagram on the screen.
</p>

<h4>Example</h4>
<p>
The example <a href=\"modelica://Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros\">
Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros</a>
is defined as
</p>
<pre>
  Plot.polesAndZeros(ss = Modelica_LinearSystems2.StateSpace(
    A=[-3, 2,-3,  4, 5,6;
        0, 6, 7,  8, 9,4;
        0, 2, 3,  0,78,6;
        0, 1, 2,  2, 3,3;
        0,13,34,  0, 0,1;
        0, 0, 0,-17, 0,0],
    B=[1,0;
       0,1;
       1,0;
       0,1;
       1,0;
       0,1],
    C=[0,0,1,0,1,0;
       0,1,0,0,1,1],
    D=[0,0;
       0,0]));
</pre>

<p>
and results in
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/polesAndZerosSS.png\"/>
</blockquote>
</html>"));
    end polesAndZeros;

    encapsulated function bodeSISO
      "Plot bode plot of the corresponding transfer function"
      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Internal;

      input StateSpace ss "State space system";
      input Integer iu=1 "Index of input";
      input Integer iy=1 "Index of output";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min=0.1
        "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max=10
        "Maximum frequency value, if autoRange = false";

      input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
      input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

      input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));
      input Boolean onFile=false
        "= true, if frequency response is stored on file as matrix [f,a,phi]" annotation(choices(checkBox=true));
      input String fileName="frequencyResponse.mat"
        "If onFile=true, file on which the frequency response will be stored"  annotation(Dialog(enable=onFile));
      input String matrixName=if Hz and not dB then "fHz_a_phiDeg" elseif
                                 Hz and dB then "fHz_adB_phiDeg" elseif
                                 not Hz and dB then "f_adB_phiDeg" else "f_a_phiDeg"
        "If onFile=true, Name of matrix on file" annotation(Dialog(enable=onFile));

    protected
      Real A[size(ss.A, 1), size(ss.A, 2)];
      Real B[size(ss.B, 1), 1];
      Real C[1, size(ss.C, 2)];
      Real D[1, 1];

      Real Eig[size(ss.A,1), 2];
      Real InvZeros[:,2];
      Real f[nPoints];
      Real a[nPoints] "Absolute value (magnitude)";
      Real phi[nPoints];
      Real gain;

      Real fap[nPoints,if onFile then 3 else 0];
      Boolean success;
    algorithm
      // Check that system has inputs and outputs
      if size(ss.B, 2) == 0 then
        print(
          "\n... Not possible to plot transfer function because system has no inputs."
           + "\n... Call of Plot.bodeSISO is ignored.\n");
        return;
      elseif size(ss.C, 1) == 0 then
        print(
          "\n... Not possible to plot transfer function because system has no outputs."
           + "\n... Call of Plot.bodeSISO is ignored.\n");
        return;
      end if;

      // Extract desired SISO system
      assert(iu <= size(ss.B, 2) and iu > 0, "Index for input is " + String(iu)
         + " which is not in [1, " + String(size(ss.B, 2)) + "].");
      assert(iy <= size(ss.C, 1) and iy > 0, "Index for output is " + String(iy)
         + " which is not in [1, " + String(size(ss.C, 1)) + "].");

      // Extract SISO system and balance it
      (,A,B,C) := Internal.balanceABC(
                    A=ss.A,
                    B=matrix(ss.B[:, iu]),
                    C=transpose(matrix(ss.C[iy, :])));
      D := matrix(ss.D[iy, iu]);

      // Compute eigenvalues and invariant zeros (as Real Matrices)
      Eig := Math.Matrices.eigenValuesAsRealMatrix(A,balance=false);
      InvZeros := StateSpace.Internal.invariantZerosWithRealMatrix(A,B,C,D);
      gain := Internal.frequencyResponseGain(A,B,C,D,InvZeros,Eig);

      // Compute frequency response values
      (f,a,phi) :=Modelica_LinearSystems2.Internal.frequencyResponse(gain,
                     InvZeros, Eig, nPoints, autoRange, f_min, f_max,
                     Hz, dB, defaultDiagram.logX);

      // Bode plot
      Internal.frequencyResponsePlot(f,a,phi,autoRange,f_min,f_max,magnitude,
                                     phase,Hz,dB,defaultDiagram,device);

      // Store frequency response values on file
      if onFile then
        fap :=[f,a,phi];
        Modelica.Utilities.Files.removeFile(fileName);
        success:=Modelica.Utilities.Streams.writeRealMatrix(
          fileName,matrixName,fap,append=false);
        if success then
          print("... Frequency response stored on file \"" +
            Modelica.Utilities.Files.fullPathName(fileName) + "\"");
        end if;
      end if;
      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>bodeSISO</b>(ss)
   or
StateSpace.Plot.<b>bodeSISO</b>(
  ss,
  iu,
  iy,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the bode-diagram of a transfer function corresponding
to the behavior of the state space system from iu'th element of the input
vector <b>u</b> to the iy'th element of the output vector <b>y</b>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0],
    B=[0.0,1.0; 1.0,1.0; -1.0,0.0],
    C=[0.0,1.0,1.0; 1.0,1.0,1.0],
    D=[1.0,0.0; 0.0,1.0])

  Integer iu=1;
  Integer iy=1;

<b>algorithm</b>
   Modelica_LinearSystems2.StateSpace.Plot.plotBodeSISO(ss, iu, iy)
//  gives:
</pre></blockquote>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</p>
</html>"));
    end bodeSISO;

    encapsulated function bodeMIMO
      "Plot bode plot of all transfer functions, corresponding to the state space system"
      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Plot;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.Internal;

      input StateSpace ss "State space system";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange[:, :] = fill(
          true,
          size(ss.C, 1),
          size(ss.B, 2)) "True, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min[:, :] = fill(
          0.1,
          size(ss.C, 1),
          size(ss.B, 2)) "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max[:, :] = fill(
          10,
          size(ss.C, 1),
          size(ss.B, 2)) "Maximum frequency value, if autoRange = false";

      input Boolean magnitude = true "= true, to plot magnitude" annotation(choices(checkBox=true));
      input Boolean phase = true "= true, to plot phase" annotation(choices(checkBox=true));

      input Real tol = 1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

      input Boolean Hz = true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB = false
        "= true, to plot magnitude in [-], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));

      input Boolean onFile = false
        "= true, if frequency response is stored on file as matrix [f,a,phi]" annotation(choices(checkBox=true));
      input String fileName = "frequencyResponse.mat"
        "If onFile=true, file on which the frequency response will be stored" annotation(Dialog(enable=onFile));
      input String matrixName = if Hz and not dB then "fHz_a_phiDeg" elseif
                                   Hz and dB then "fHz_adB_phiDeg" elseif
                                   not Hz and dB then "f_adB_phiDeg" else "f_a_phiDeg"
        "If onFile=true, prefix name of matrix on file" annotation(Dialog(enable=onFile));

    protected
      Real A[    size(ss.A, 1), size(ss.A, 2)];
      Real Bfull[size(ss.B, 1), size(ss.B, 2)];
      Real Cfull[size(ss.C, 1), size(ss.C, 2)];

      Real B[size(ss.B, 1), 1];
      Real C[1, size(ss.C, 2)];
      Real D[1,1];

      Real Eig[size(ss.A,1), 2];
      Real InvZeros[:,2];
      Real f[nPoints] "Frequency";
      Real a[nPoints] "Absolute value (magnitude)";
      Real phi[nPoints] "Phase";
      Real gain;

      Real fap[nPoints,if onFile then 3 else 0];
      Boolean success;

      String yNames[size(ss.C, 1)];
      String uNames[size(ss.B, 2)];

      Plot.Records.Diagram diagram2 = defaultDiagram;
    algorithm
      // Check that system has inputs and outputs
      if size(ss.B, 2) == 0 then
        print("\n... Not possible to plot transfer function because system has no inputs."
           + "\n... Call of Plot.bodeMIMO is ignored.\n");
        return;
      elseif size(ss.C, 1) == 0 then
        print("\n... Not possible to plot transfer function because system has no outputs."
           + "\n... Call of Plot.bodeMIMO is ignored.\n");
        return;
      end if;

      // generate headings
      for i1 in 1:size(ss.B, 2) loop
        uNames[i1] := if ss.uNames[i1] == "" then "u" + String(i1) else ss.uNames[i1];
      end for;
      for i1 in 1:size(ss.C, 1) loop
        yNames[i1] := if ss.yNames[i1] == "" then "y" + String(i1) else ss.yNames[i1];
      end for;

      // Balance system
      (,A,Bfull,Cfull) := Internal.balanceABC(A=ss.A, B=ss.B, C=ss.C);

      // Compute eigen values (matrix A balanced already)
      Eig := Math.Matrices.eigenValuesAsRealMatrix(A,balance=false);

      // Remove output file, if onFile=true
      if onFile then
         Modelica.Utilities.Files.removeFile(fileName);
      end if;

      // Perform computation from every input to every output
      for i1 in 1:size(ss.C, 1) loop
        for i2 in 1:size(ss.B, 2) loop
          // Compute zeros
          B := matrix(Bfull[:, i2]);
          C := transpose(matrix(Cfull[i1, :]));
          D := matrix(ss.D[i1, i2]);
          InvZeros := StateSpace.Internal.invariantZerosWithRealMatrix(A,B,C,D);
          gain := Internal.frequencyResponseGain(A,B,C,D,InvZeros,Eig);

          // Compute frequency response values
          (f,a,phi) := Modelica_LinearSystems2.Internal.frequencyResponse(
                         gain, InvZeros, Eig, nPoints, autoRange[i1, i2],
                         f_min[i1, i2], f_max[i1, i2], Hz, dB, defaultDiagram.logX);

          // Bode plot
          diagram2.heading := defaultDiagram.heading + "  " + uNames[i2] + " -> " + yNames[i1];
          Internal.frequencyResponsePlot(f,a,phi,autoRange[i1, i2],
                                         f_min[i1, i2], f_max[i1, i2],magnitude,
                                         phase,Hz,dB,diagram=diagram2,
                                         device=device);

          // Store result optionally on file
          if onFile then
            fap := [f,a,phi];
            success := Modelica.Utilities.Streams.writeRealMatrix(
              fileName, matrixName + "_" + uNames[i2] + "_" + yNames[i1], fap, append=true);
          end if;
        end for;
      end for;

      if onFile and success then
        print("... Frequency response stored on file \"" +
          Modelica.Utilities.Files.fullPathName(fileName) + "\"");
      end if;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>bodeMIMO</b>(ss)
   or
StateSpace.Plot.<b>bodeMIMO</b>(
  ss,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude,
  phase,
  tol,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>(),
  Hz,
  dB,
  onFile,
  fileName,
  matrixName)
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0],
    B=[0.0,1.0; 1.0,1.0; -1.0,0.0],
    C=[0.0,1.0,1.0],
    D=[1.0,0.0])

<b>algorithm</b>
   Modelica_LinearSystems2.StateSpace.Plot.plotBodeMIMO(ss)
//  gives:
</pre></blockquote>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude2.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase2.png\">
</p>
</html>"));
    end bodeMIMO;

    encapsulated function timeResponse
      "Plot the time response of the system. The response type is selectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Modelica.SIunits.Time dt=0 "Sample time";
      input Modelica.SIunits.Time tSpan=0 "Simulation time span";

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step;

      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
        annotation (Dialog, choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Time response"));

    protected
      Plot.Records.Curve curve;
      Integer i1;
      Integer i2;
      Plot.Records.Diagram diagram2[size(ss.C, 1)];

      Real y[:, size(ss.C, 1), if response == TimeResponse.Initial then 1 else
        size(ss.B, 2)]
        "Output response: (number of samples) x (number of outputs) x (number of inputs)";
      Real t[:] "Time vector: (number of samples)";
      Real x[:, size(ss.A, 1), if response == TimeResponse.Initial then 1 else
        size(ss.B, 2)]
        "State trajectories: (number of samples) x (number of states) x (number of inputs)";
      String yNames[size(ss.C, 1)];
      String uNames[size(ss.B, 2)];
      Integer loops=if response == TimeResponse.Initial then 1 else size(ss.B,
          2);

    algorithm
      (y,t,x) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
            sc=ss,
            dt=dt,
            tSpan=tSpan,
            response=response,
            x0=x0);

      // generate headings
      for i1 in 1:size(ss.B, 2) loop
        uNames[i1] := if ss.uNames[i1] == "" then "u" + String(i1) else ss.uNames[
          i1];
      end for;
      for i1 in 1:size(ss.C, 1) loop
        yNames[i1] := if ss.yNames[i1] == "" then "y" + String(i1) else ss.yNames[
          i1];
      end for;

      for i2 in 1:loops loop
        for i1 in 1:size(ss.C, 1) loop
          curve := Plot.Records.Curve(
                x=t,
                y=y[:, i1, i2],
                autoLine=true);

          diagram2[i1] := defaultDiagram;
          diagram2[i1].curve := {curve};
          diagram2[i1].heading := if response == TimeResponse.Initial then
            defaultDiagram.heading + " " + yNames[i1] else defaultDiagram.heading
             + "  " + uNames[i2] + " -> " + yNames[i1];
          diagram2[i1].yLabel := yNames[i1];

        end for;

        if subPlots then
          Plot.diagramVector(diagram2, device);
        else
          for i1 in 1:size(ss.C, 1) loop
            Plot.diagram(diagram2[i1], device);
          end for;
        end if;
      end for;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>timeResponse</b>(ss);
   or
StateSpace.Plot.<b>timeResponse</b>(
  ss,
  dt,
  tSpan,
  response,
  x0,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the time response of a state space system. The character of the time response if defined by the input <tt>response</tt>, i.e. Impulse, Step, Ramp or Initial.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
    B=[1.0; 1.0; 0.0],
    C=[0.0,1.0,1.0],
    D=[0.0])

  Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Plot.timeResponse(ss, response=response)
// gives:
</pre></blockquote>
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/timeResponseSS.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.initialResponse\">initial</a>
</p>
</html>"));
    end timeResponse;

    encapsulated function impulse "Impulse response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Modelica.SIunits.Time dt=0 "Sample time";
      input Modelica.SIunits.Time tSpan=0 "Simulation time span";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
        annotation (Dialog, choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Impulse response"));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "type of time response";
    algorithm

      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
            ss=ss,
            dt=dt,
            tSpan=tSpan,
            response=response,
            x0=x0,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>impulse</b>(ss);
   or
StateSpace.Plot.<b>impulse</b>(
  ss,
  dt,
  tSpan,
  x0,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the impulse responses of a state space system
for each system corresponding to the transition matrix. It is based on <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
    B=[1.0; 1.0; 0.0],
    C=[0.0,1.0,1.0],
    D=[0.0])

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Plot.impulse(ss)
// gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/impulseResponseSS.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
    end impulse;

    encapsulated function step "Step response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Modelica.SIunits.Time dt=0 "Sample time";
      input Modelica.SIunits.Time tSpan=0 "Simulation time span";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
        annotation (Dialog, choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Step response"));

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";

    algorithm
      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
            ss=ss,
            dt=dt,
            tSpan=tSpan,
            response=response,
            x0=x0,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>step</b>(ss);
   or
StateSpace.Plot.<b>step</b>(
  ss,
  dt,
  tSpan,
  x0,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the step responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
    B=[1.0; 1.0; 0.0],
    C=[0.0,1.0,1.0],
    D=[0.0])

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Plot.step(ss, tSpan=3)
// gives:
</pre></blockquote>
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/stepResponseSS.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
    end step;

    encapsulated function ramp "Ramp response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Modelica.SIunits.Time dt=0 "Sample time";
      input Modelica.SIunits.Time tSpan=0 "Simulation time span";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
        annotation (Dialog, choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Ramp response"));

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";

    algorithm
      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
            ss=ss,
            dt=dt,
            tSpan=tSpan,
            response=response,
            x0=x0,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>ramp</b>(ss);
   or
StateSpace.Plot.<b>ramp</b>(
  ss,
  dt,
  tSpan,
  x0,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the ramp responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
    B=[1.0; 1.0; 0.0],
    C=[1.0,1.0,1.0],
    D=[0.0])

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Plot.ramp(ss)
// gives:
</pre></blockquote>
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/rampResponseSS.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
    end ramp;

    encapsulated function initialResponse "Initial condition response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Modelica.SIunits.Time dt=0 "Sample time";
      input Modelica.SIunits.Time tSpan=0 "Simulation time span";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
        annotation (Dialog, choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading=
             "Initial response"));

      input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";
    algorithm

      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
            ss=ss,
            dt=dt,
            tSpan=tSpan,
            response=response,
            x0=x0,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>initial</b>(ss);
   or
StateSpace.Plot.<b>initial</b>(
  ss,
  dt,
  tSpan,
  x0,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the initial responses of a state space system for the initial state vector x0 for each system corresponding to the transition matrix. It is based on <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
    B=[1.0; 1.0; 0.0],
    C=[0.0,1.0,1.0],
    D=[0.0])

  Real x0={1,0.5,0.5};

<b>algorithm</b>
  Modelica_LinearSystems2.StateSpace.Plot.initial(ss, x0=x0)
// gives:
</pre></blockquote>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/initialResponseSS.png\">
</blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>,
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>
</p>
</html>"));
    end initialResponse;

  end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of StateSpace data record"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function toZerosAndPoles
      "Generate a zeros-and-poles representation from a SISO state space representation"

      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      //protected
      //  input Boolean cancel=true "false to hinder cancellation";// is not fully realized
    public
      output ZerosAndPoles zp;

    protected
      StateSpace ssm=if size(ss.A, 1) > 0 then
          StateSpace.Transformation.toIrreducibleForm(ss, tol) else
          StateSpace(ss.D[1, 1]);
      Complex Poles[:];
      Complex Zeros[:];

      Real gain;

      Complex frequency;
      Complex Gs;
      Real As[:, :];
      Real pk;
      Integer i;
      Integer k;
      Boolean h;
      Real v;

      function getReOutsidePolesZeros
        "Get p on the real axis so that there is a minimum distance to all poles and zeros"
        input Complex Poles[:];
        input Complex Zeros[:];
        input Real minimumDistance=0.1;
        output Real p
          "Value on real axis > 0.0, so that poles.re and zeros.re have a minimumDistance to it";
        /* Most systems have no or only a few unstable poles or zeros.
    Searching for a suitable p is therefore fastest when searching 
    only in the unstable region, that is p > 0.0
    */
      protected
        Integer nVec=size(Poles, 1) + size(Zeros, 1);
        Real vec[:];
        Real vecSorted[:];
        Integer i;
        Integer iMax;
        Real small = minimumDistance*1e-6;
      algorithm
        i := 0;
        vec := zeros(nVec);
        for j in 1:size(Poles, 1) loop
          if Poles[j].re > small then
            i := i + 1;
            vec[i] := Poles[j].re;
          end if;
        end for;
        for j in 1:size(Zeros, 1) loop
          if Zeros[j].re > small then
            i := i + 1;
            vec[i] := Zeros[j].re;
          end if;
        end for;
        iMax := i;
        if iMax == 0 then
          p := minimumDistance;
          return;
        end if;

        vec := vec[1:iMax];
        vecSorted := Modelica.Math.Vectors.sort(vec);

        // Find p, so that vecSorted[i+1] - vecSorted[i] > 2*minimumDistance
        if vecSorted[1] >= 2*minimumDistance then
           p := minimumDistance;
           return;
        end if;

        i :=1;
        while i <= iMax loop
          if i == iMax then
            p := vecSorted[i] + minimumDistance;
            break;
          elseif vecSorted[i + 1] - vecSorted[i] > 2*minimumDistance then
            p := vecSorted[i] + minimumDistance;
            break;
          else
            i := i + 1;
          end if;
        end while;
      end getReOutsidePolesZeros;

    algorithm
      if Modelica.Math.Vectors.length(ssm.B[:, 1]) > 0 and
          Modelica.Math.Vectors.length(ssm.C[1, :]) > 0 then
        Poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.eigenValues_dhseqr(ssm.A);
        //ssm.A is of upper Hessenberg form
        Zeros := StateSpace.Analysis.invariantZeros(ssm);

        if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
          assert(size(ss.B, 2) == 1,
            " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
             + String(size(ss.B, 2)) + " instead of 1");
          assert(size(ss.C, 1) == 1,
            " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
             + String(size(ss.C, 1)) + " instead of 1");
        end if;
        zp := ZerosAndPoles(z=Zeros, p=Poles, k=1);

        v := getReOutsidePolesZeros(Poles, Zeros);
        frequency := Complex(v);
        Gs := ZerosAndPoles.Analysis.evaluate(zp, frequency);
        As := -ssm.A;
        for i in 1:size(As, 1) loop
          As[i, i] := As[i, i] + frequency.re;
        end for;

        pk := StateSpace.Internal.partialGain(As, ssm.B[:, 1]);
        gain := (ssm.C[1, size(As, 1)]*pk + ss.D[1, 1])/Gs.re;
        zp := ZerosAndPoles(
              z=Zeros,
              p=Poles,
              k=gain);
      else
        zp := ZerosAndPoles(
              z=fill(Complex(0), 0),
              p=fill(Complex(0), 0),
              k=scalar(ss.D));

      end if;
      zp.uName := ss.uNames[1];
      zp.yName := ss.yNames[1];

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = StateSpace.Conversion.<b>toZerosAndPoles</b>(ss)
</pre> </blockquote>

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
of a system from state space representation using the transformation algorithm described in [1].
The uncontrollable and unobservable parts are isolated and the eigenvalues and invariant zeros of the controllable and observable sub system are calculated.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-1.0, 0.0, 0.0;
          0.0,-2.0, 0.0;
          0.0, 0.0,-3.0],
    B = [1.0;
         1.0;
         0.0],
    C = [1.0,1.0,1.0],
    D = [0.0]);

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);
//                s + 1.5
//   zp = 2 -----------------
             (s + 1)*(s + 2)
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2011-07-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Improved frequency calculation</td>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
    end toZerosAndPoles;

    function toTransferFunction
      "Generate a transfer function from a SISO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      output TransferFunction tf;

    protected
      ZerosAndPoles zp;

    algorithm
      zp := toZerosAndPoles(ss, tol);
      tf := ZerosAndPoles.Conversion.toTransferFunction(zp);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tf = StateSpace.Conversion.<b>toTransferFunction</b>(ss)
</pre> </blockquote>

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
The algorithm uses <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles\">StateSpace.Conversion.toZerosAndPoles</a> to convert the state space system into a zeros and poles representation first and after that <a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction\">ZerosAndPoles.Conversion.toTransferFunction</a> to generate the transfer function.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-1.0, 0.0, 0.0;
          0.0,-2.0, 0.0;
          0.0, 0.0,-3.0],
    B = [1.0;
         1.0;
         0.0],
    C = [1.0,1.0,1.0],
    D = [0.0]);

<b>algorithm</b>
  tf:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);
//             2*s + 3
//   tf =  -----------------
             s^2 + 3*s + 2
</pre></blockquote>
</html>", revisions="<html>
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
    end toTransferFunction;

    encapsulated function toZerosAndPolesMIMO
      "Generate a zeros-and-poles representation from a MIMO state space representation"
      import Modelica;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      output ZerosAndPoles zp[size(ss.C, 1), size(ss.B, 2)];

    protected
      StateSpace ss_siso(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), 1],
        redeclare Real C[1, size(ss.C, 2)],
        redeclare Real D[1, 1]);

      Integer ny=size(ss.C, 1);
      Integer nu=size(ss.B, 2);

    algorithm
      for ic in 1:ny loop
        for ib in 1:nu loop
          ss_siso := StateSpace(
                A=ss.A,
                B=matrix(ss.B[:, ib]),
                C=transpose(matrix(ss.C[ic, :])),
                D=matrix(ss.D[ic, ib]));
          zp[ic, ib] := StateSpace.Conversion.toZerosAndPoles(ss_siso, tol);
        end for;
      end for;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = StateSpace.Conversion.<b>toZerosAndPolesMIMO</b>(ss)
</pre> </blockquote>

<h4>Description</h4>
<p>
Computes a matrix of ZerosAndPoles records
</p>
<blockquote><pre>
          product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
zp = k * ---------------------------------------------------------
          product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
<p>
of a system from state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each input-output pair.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-1.0, 0.0, 0.0;
          0.0,-2.0, 0.0;
          0.0, 0.0,-3.0],
    B = [0.0, 1.0;
         1.0, 1.0;
        -1.0, 0.0],
    C = [0.0, 1.0, 1.0;
         1.0, 1.0, 1.0],
    D = [1.0, 0.0;
         0.0, 1.0]);

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);

// zp = [(s^2 + 5*s + 7)/( (s + 2)*(s + 3) ), 1/(s + 2);
         1/( (s + 2)*(s + 3) ), 1*(s + 1.38197)*(s + 3.61803)/( (s + 1)*(s + 2) )]
</pre></blockquote>
<p>
i.e.
</p>
<blockquote><pre>
         |                                                   |
         |    (s^2+5*s+7)                    1               |
         | -----------------               -----             |
         |  (s + 2)*(s + 3)                (s+2)             |
  tf  =  |                                                   |
         |        1             (s + 1.38197)*(s + 3.61803)  |
         | -------------       ----------------------------- |
         | (s + 2)*(s + 3)            (s + 1)*(s + 2)        |
         |                                                   |
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
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
    end toZerosAndPolesMIMO;

    function toTransferFunctionMIMO
      "Generate a transfer function of a MIMO system from state space representation"

      import Modelica;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      output TransferFunction tf[size(ss.C, 1), size(ss.B, 2)]
        "Matrix of transfer function objects";

    protected
      ZerosAndPoles zp[:, :];
      parameter Integer m=size(ss.B, 2);
      parameter Integer p=size(ss.C, 1);

    algorithm
      zp := StateSpace.Conversion.toZerosAndPolesMIMO(ss, tol);
      for i1 in 1:m loop
        for i2 in 1:p loop
          tf[i2, i1] := ZerosAndPoles.Conversion.toTransferFunction(zp[i2, i1]);
        end for;
      end for;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tf = StateSpace.Conversion.<b>toTransferFunctionMIMO</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a matrix of TransferFunction records
</p>
<blockquote><pre>
        n_i(s)     b0_i + b1_i*s + ... + bn_i*s^n
tf_i = -------- = --------------------------------
        d_i(s)     a0_i + a1_i*s + ... + an_i*s^n
</pre></blockquote>
<p>
with repetitive application of <a href=\"modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">Conversion.toTransferFunction</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-1.0, 0.0, 0.0;
          0.0,-2.0, 0.0;
          0.0, 0.0,-3.0],
    B = [0.0, 1.0;
         1.0, 1.0;
        -1.0, 0.0],
    C = [0.0, 1.0, 1.0;
         1.0, 1.0, 1.0],
    D = [1.0, 0.0;
         0.0, 1.0]);

<b>algorithm</b>
  zp:=Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles(ss);

// zp = [(s^2 + 5*s + 7)/(s^2 + 5*s + 6), 1/(s + 2);
         1/(s^2 + 5*s + 6), (1*s^2 + 5*s + 5)/(s^2 + 3*s + 2)]
</pre></blockquote>
<p>
i.e.
</p>
<blockquote><pre>
         |                                                   |
         |    (s^2+5*s+7)                    1               |
         | -----------------               -----             |
         |  (s + 2)*(s + 3)                (s+2)             |
  tf  =  |                                                   |
         |        1             (s + 1.38197)*(s + 3.61803)  |
         | -------------       ----------------------------- |
         | (s + 2)*(s + 3)            (s + 1)*(s + 2)        |
         |                                                   |
</pre></blockquote>
</html>", revisions="<html>
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
    end toTransferFunctionMIMO;

  end Conversion;

  encapsulated package Transformation
    "Package of functions for state space similarity transformations"
    import Modelica;
    extends Modelica.Icons.Package;

    encapsulated function toSimilarForm
      "Perform the similarity transformation z = Tx (or x = inv(T)z) which leads to Az=T*A*inv(T), Bz=T*B, Cz=C*inv(T), Dz=D (or Az=inv(T)*A*T, Bz=inv(T)B, Cz=C*T, Dz=D)"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      input Real T[size(ss.A, 2), size(ss.A, 1)]=identity(size(ss.A, 1))
        "Transformation matrix";
      input Boolean inverted=false
        "Is false (default) for transformation z = Tx, true for x = Tz"
        annotation (choices(checkBox=true));

      output StateSpace tss(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);

    algorithm
      if inverted then
        tss.A := Matrices.solve2(T, ss.A*T);
        tss.B := Matrices.solve2(T, ss.B);
        tss.C := ss.C*T;
        tss.D := ss.D;
      else
        tss.A := transpose(LAPACK.dgesvx(T, transpose(T*ss.A)));
        tss.B := T*ss.B;
        tss.C := transpose(LAPACK.dgesvx(T, transpose(ss.C)));
        tss.D := ss.D;
      end if;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toSimilarForm</b>(ss, T, inverted)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates a similar state space system, i.e.
</p>
<blockquote><pre>
der(z) = T*A*inv(T)*z + T*B*u
     y = C*T*z + D*u
</pre></blockquote>
<p>
if inverted==false and
</p>
<blockquote><pre>
der(z) = inv(T)*A*T*z + inv(T)*B*u
     y = C*inv(T)*z + D*u
</pre></blockquote>
<p>
if inverted=true. Matrix T has to be invertible. The transformed system has the same eigenvalues.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
     A=[-1, 1; 0, -2],
     B=[1; 0],
     C=[0, 1],
     D=[0]);

  Real T[2,2]=[1, 1;0, sqrt(2)];

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm(ss, T, false);
//  tss=StateSpace(
      A=[-1, 0; 0, -2],
      B=[1; 0],
      C=[0, 1/sqrt82)],
      D=[0])
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Analysis.analysis\">State space analysis</a>
</p>
</html>"));
    end toSimilarForm;

    encapsulated function toObservabilityForm
      "Perform the similarity transformation to the obervabillity canonical form"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      output StateSpace tss(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);

    protected
      Real V[size(ss.A, 1), size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

      Integer nx=size(ss.A, 1);

    algorithm
      assert(size(ss.C, 1) == 1 and size(ss.B, 2) == 1,
        "Calculation of controllable form fails for systems with more than 1 inputs or outputs");
      assert(Modelica_LinearSystems2.StateSpace.Analysis.isObservable(ss),
        "transformation ist not realizable since the system ist not obersvable");

      V[:, 1] := Modelica.Math.Matrices.solve(
        StateSpace.Analysis.observabilityMatrix(ss), vector([fill(
            0,
            1,
            nx - 1), 1]));

      for i in 2:nx loop
        V[:, i] := ss.A*V[:, i - 1];
      end for;

      tss := StateSpace.Transformation.toSimilarForm(
            ss,
            V,
            inverted=true);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toObservabilityForm</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the observability form of a SISO state space system, i.e.
</p>
<blockquote><pre>
tss:
der(z) = inv(T)*A*T*z + inv(T)*B*u
     y = C*inv(T)*z + D*u
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
T = [C; C*A; ...; C*A^(n-1)]
</pre></blockquote>
<p>
is the observability matrix of the original state space system.
In comparison to the corresponding transfer function
</p>
<blockquote><pre>
        b0 + b1*s + ... + bn*s^n
G(s) = --------------------------
        a0 + a1*s + ... + an*s^n
</pre></blockquote>
<p>
the canonical observability form is
</p>
<blockquote><pre>
    | 0   0   ...   0   -a0   |        | b0   - a0*bn   |
    | 1   0   ...   0   -a1   |        | b1   - a1*bn   |
A = | 0   1   ...   0   -a2   |,   B = |     ...        |
    |... ...  ...  ...  -a3   |        | bn-2 - an-2*bn |
    | 0  ...  ...   1   -an-1 |        | bn-1 - an-1*bn |

C = [0, 0, ..., 1],                D = [bn]
</pre></blockquote>
<p>
Matrix <code>T</code> has to be invertible, i.e. the system has to be observable.
The transformed system has the same eigenvalues.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1; 1, -2],
    B=[1; 0],
    C=[1, 1],
    D=[2]);

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toObservabilityForm(ss);
//  tss=StateSpace(
      A=[0, -1; 1, -3],
      B=[3; 1],
      C=[0, 1],
      D=[2])
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toControllabilityForm\">toControllabilityForm</a>
</p>
</html>", revisions="<html>
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
    end toObservabilityForm;

    encapsulated function toControllabilityForm
      "Perform the similarity transformation to the controllability canonical form"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      output StateSpace tss(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);

    protected
      Real V[size(ss.A, 1), size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

      Integer nx=size(ss.A, 1);

    algorithm
      assert(size(ss.C, 1) == 1 and size(ss.B, 2) == 1,
        "Calculation of controllable form fails for systems with more than 1 inputs or outputs");
      assert(Modelica_LinearSystems2.StateSpace.Analysis.isControllable(ss),
        "transformation ist not realizable since the system ist not controllable");

      V[1, :] := Modelica.Math.Matrices.solve(transpose(
        StateSpace.Analysis.controllabilityMatrix(ss)), vector([fill(
            0,
            nx - 1,
            1); 1]));

      for i in 2:nx loop
        V[i, :] := V[i - 1, :]*ss.A;
      end for;

      tss := StateSpace.Transformation.toSimilarForm(ss, V);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toControllabilityForm</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the controllability form of a SISO state space system, i.e.
</p>
<blockquote><pre>
tss:
der(z) = T*A*inv(T)*z + T*B*u
     y = C*T*z + D*u
</pre></blockquote>
<p>
with
</p>
<blockquote><pre>
T = [B, A*B,..., A^(n-1)*B]
</pre></blockquote>
<p>
is the observability matrix of the original state space system.
In comparison to the corresponding transfer function
</p>
<blockquote><pre>
        b0 + b1*s + ... + bn*s^n
G(s) = --------------------------
        a0 + a1*s + ... + an*s^n
</pre></blockquote>
<p>
the canonical observability form is
</p>
<blockquote><pre>
    | 0   1   0   ...   0     0   |                        | 0 |
    |  0     0     0    0     0   |                        | 0 |
A = | ...   ...   ...  ...   ...  |,                   B = |...|
    |  0     0     0    0     0   |                        | 0 |
    | -a0   -a1   -a2  ...  -an-1 |                        | 1 |

C = [ b0 - bn*a0, b1 - bn*a1, ..., bn-1 - bn*an-1],    D = [bn]
</pre></blockquote>
<p>
Matrix T has to be invertible, i.e. the system has to be controllable. The transformed system has the same eigenvalues.
</p>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1, 1; 1, -2],
      B=[1; 0],
      C=[1, 1],
      D=[2]);

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toControllabilityForm(ss);
//  tss=StateSpace(
      A=[0, 1; -1, -3],
      B=[0; 1],
      C=[3, 1],
      D=[2])
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toObservabilityForm\">toObservabilityForm</a>
</p>
</html>", revisions="<html>
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
    end toControllabilityForm;

    encapsulated function toDiagonalForm
      "Perform the similarity transformation with the (real) inverse right eigenvector matrix of the system, that lead to the Jordan canonical form for single eigenvalues"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      output StateSpace tss(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);

    protected
      Real V[size(ss.A, 1), size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

    algorithm
      (,,,V,) := LAPACK.dgeev(ss.A);

      tss := StateSpace.Transformation.toSimilarForm(
            ss,
            V,
            inverted=true);

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toDiagonalForm</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the diagonal form of a SISO state space system, i.e.
</p>
<blockquote><pre>
tss:
der(z) = inv(T)*A*T*z + inv(T)*B*u
     y = C*inv(T)*z + D*u
</pre></blockquote>
<p>
Matrix T has to be diagonalizable, i.e. the algebraic and geometric multiplicities of an eigenvalue must coincide. The diagonal entries of the new system matrix tss.<b>A</b> are the eigenvalues off the systemmatrix ss.<b>A</b>.
</p>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1, 1; 0, -2],
      B=[1; 0],
      C=[1, 1],
      D=[2]);

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toDiagonalForm(ss);
//  tss=StateSpace(
      A=[-1, 0; 0, -2],
      B=[1; 0],
      C=[1, 0],
      D=[0])
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>
</p>
</html>"));
    end toDiagonalForm;

    encapsulated function toBalancedForm
      "Perform the similarity transformation to a balanced form (to make further numerical computations on the StateSpace system more reliable)"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output StateSpace ssBalanced(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]) "Balanced ss";
    algorithm
      (,ssBalanced.A, ssBalanced.B, ssBalanced.C) :=
        Modelica_LinearSystems2.Internal.balanceABC(
            ss.A,
            ss.B,
            ss.C);
      ssBalanced.D :=ss.D;
      ssBalanced.uNames := ss.uNames;
      ssBalanced.yNames := ss.yNames;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ssBalanced = StateSpace.Transformation.<b>toBalancedForm</b>(ss);
</pre></blockquote>

<h4>Description</h4>

<p>
Balancing a linear dynamic system in state space form ss means to find a 
state transformation x_new = T*x = diagonal(scale)*x
so that the transformed system is better suited for numerical algorithms.
In more detail:
</p>

<p>
This function performs a similarity transformation with T=diagonal(scale) such that S_scale
</p>

<pre>             |inv(T)*ss.A*T, inv(T)*ss.B|
   S_scale = |                          |
             |       ss.C*T,     0      |
</pre>

<p>
has a better condition as system matrix S
</p>

<pre>       |ss.A, ss.B|
   S = |          |
       |ss.C, 0   |
</pre>
that is, conditionNumber(S_scale) &le; conditionNumber(S). The elements of vector scale
are multiples of 2 which means that this function does not introduce round-off errors.
</p>


<h4>Example</h4>

<blockquote>
<pre>import Modelica.Math.Matrices.norm;
ss = Modelica_LinearSystems2.StateSpace(A=[1, -10,  1000; 0.01,  0,  10; 0.005,  -0.01,  10], 
                                        B=[100, 10; 1,0; -0.003, 1], 
                                        C=[-0.5, 1, 100], 
                                        D=[0,0]);
sb = Modelica_LinearSystems2.StateSpace.Transformation.toBalancedForm(ss);

-> Results in: 
norm(ss.A) = 1000.15, norm(ss.B) = 100.504, norm(ss.C) = 100.006
norm(sb.A) = 10.8738, norm(sb.B) = 16.0136, norm(sb.C) = 10.2011
</pre>
</blockquote>

<p>
The algorithm is taken from
</p>
<dl>
<dt>H. D. Joos, G. Gr&uuml;bel:
<dd><b>RASP'91 Regulator Analysis and Synthesis Programs</b><br>
    DLR - Control Systems Group 1991
</dl>
<p>
which is based on the <code>balance</code> function from EISPACK.
</p>
</html>"));
    end toBalancedForm;

    encapsulated function toIrreducibleForm
      "Calculate a minimal controllable and observable block Hessenberg realization of a given SISO state-space representation "

      // test of SISO has to be added
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssm1=
          StateSpace.Internal.reducedCtrSystem(ss, tol);
      Integer nx=size(ss.A, 1);
      Integer rankQ=ssm1.r;
      StateSpace ss2=StateSpace(
              A=transpose(ssm1.A[nx - rankQ + 1:nx, nx - rankQ + 1:nx]),
              B=transpose(ssm1.C[:, nx - rankQ + 1:nx]),
              C=transpose(ssm1.B[nx - rankQ + 1:nx, :]),
              D=ssm1.D);
      Integer nx2=ssm1.r;
      Modelica_LinearSystems2.Internal.StateSpaceR ssm2=
          StateSpace.Internal.reducedCtrSystem(ss2, tol);
      Integer rankQ2=ssm2.r;
    public
      output StateSpace ssm3(
        redeclare Real A[rankQ2, rankQ2],
        redeclare Real B[rankQ2, size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), rankQ2],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);
    algorithm
      ssm3 := StateSpace(
            A=transpose(ssm2.A[nx2 - rankQ2 + 1:nx2, nx2 - rankQ2 + 1:nx2]),
            B=transpose(ssm2.C[:, nx2 - rankQ2 + 1:nx2]),
            C=transpose(ssm2.B[nx2 - rankQ2 + 1:nx2, :]),
            D=(ssm2.D));

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
tss = StateSpace.Transformation.<b>toIrreducibleForm</b>(ss)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates a minimal controllable and observable block Hessenberg realization for a given state-space representation.
Therefore, all uncontrollable and unobservable modes are removed by performing orthogonal similarity transformations as described in [1].
</p>
<p>
This function is called to compute transfer functions of state space representations as described in [1]. Look at [1] for further details
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A = [-4.5,  1.5,   4.0;
         -4.0,  1.0,   4.0;
         -1.5, -0.5,   1.0],
    B = [  1; 0; 1 ],
    C = [1,  0,  0],
    D = [0]);

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(ss);
//  tss=StateSpace(
      A=[-0.5],
      B=[-sqrt(0.5)],
      C=[-1/sqrt(0.5)1],
      D=[0])
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
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
    end toIrreducibleForm;

    encapsulated function toStaircaseForm
      "Transforms a state space system to upper staircase form"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD "Method for staircase algorithm";

      output StateSpace ss_sc(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)]);

    algorithm
      if method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD then
        (,ss_sc) := StateSpace.Internal.staircaseSVD(ss);
      else
        (,ss_sc) := Modelica_LinearSystems2.StateSpace.Internal.staircaseQR(ss);
      end if;

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss_sc = StateSpace.Transformation.<b>toStaircaseForm</b>(ss, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the upper staircase form state space system.
</p>

<h4>Example</h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
   A=[17.0,   24.0,    1.0,    8.0,   15.0;
      23.0,    5.0,    7.0,   14.0,   16.0;
       4.0,    6.0,   13.0,   20.0,   22.0;
      10.0,   12.0,   19.0,   21.0,    3.0;
      11.0,   18.0,   25.0,    2.0,    9.0],
   B=[-1.0,   -4.0;
       4.0,    9.0;
      -9.0,  -16.0;
      16.0,   25.0;
     -25.0,  -36.0],
   C=[1, 0, 1, 0, 0;
      0, 1, 0, 1, 1],
   D=[0, 0;
      0, 0]);

<b>algorithm</b>
  ss_sc:=Modelica_LinearSystems2.StateSpace.Transformation.toStaircaseForm(ss);
  ss_sc=StateSpace(
    A=[-1, 0; 0, -2],
    B=[1; 0],
    C=[1, 0],
    D=[0])
</pre></blockquote>
</html>", revisions="<html>
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
    end toStaircaseForm;

    encapsulated function extract
      "Extract input/output related subsystems from state space system record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Integer outputIndex[:]={1} "Vector of subsystem outputs indices";
      input Integer inputIndex[:]={1} "Vector of subsystem inputs indices";
      output StateSpace subSc(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(inputIndex, 1)],
        redeclare Real C[size(outputIndex, 1), size(ss.C, 2)],
        redeclare Real D[size(outputIndex, 1), size(inputIndex, 1)])
        "Subsystem state space record";

    algorithm
      subSc.A := ss.A;
      subSc.B := ss.B[:, inputIndex];
      subSc.C := ss.C[outputIndex, :];
      subSc.D := ss.D[outputIndex, inputIndex];

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
subsystem = StateSpace.Transformation.<b>extract</b>(ss, outputIndex, inputIndex)
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the subsystem of a state space system corresponding to the indices in outputIndex and inputIndex, i.e.
</p>
<blockquote><pre>
subsystem.A = ss.A;
subsystem.B = ss.B[:, inputIndex];
subsystem.C = ss.C[outputIndex, :];
subsystem.D = ss.D[outputIndex, inputIndex];
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1, 1, 2; 0, -2, 3;-3, 2, 1],
    B=[1, 0; 0, 1; 1, 0],
    C=[1, 1, 0; 0, 1, 1],
    D=[0, 0; 0,0]);

  Integer outputIndex={1, 2};
  Integer inputIndex={2}

<b>algorithm</b>
  tss:=Modelica_LinearSystems2.StateSpace.Transformation.extract(ss, outputIndex, inputIndex);
//  tss=StateSpace(
      A=[-1, 1, 2; 0, -2, 3;-3, 2, 1],
      B=[0; 1; 0],
      C=[1, 1, 0; 0, 1, 1],
      D=[0; 0])
</pre></blockquote>
</html>"));
    end extract;

  end Transformation;

  encapsulated package Import
    "Package of functions to generate a StateSpace data record from imported data"
    extends Modelica.Icons.Package;
    import Modelica;

    encapsulated function fromFile
      "Read a StateSpace data record from mat-file"
      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

      input String fileName = "dslin.mat" "Name of the state space system data file"
        annotation (
          Dialog(
            loadSelector(
              filter="MAT files (*.mat);; All files (*.*)",
              caption="State space system data file")));
      input String matrixName = "ABCD"
        "Name of the state space system matrix (default is \"ABCD\") in the fileName"
        annotation (Dialog);
    protected
      Integer xuy[3] = ReadSystemDimension(fileName, matrixName) annotation(__Dymola_allowForSize=true);
      Integer nx = xuy[1] annotation(__Dymola_allowForSize=true);
      Integer nu = xuy[2] annotation(__Dymola_allowForSize=true);
      Integer ny = xuy[3] annotation(__Dymola_allowForSize=true);

    public
      output StateSpace result(
        redeclare Real A[nx, nx],
        redeclare Real B[nx, nu],
        redeclare Real C[ny, nx],
        redeclare Real D[ny, nu]) "Model read from file";

    protected
      Real ABCD[nx + ny, nx + nu] = Streams.readRealMatrix(
        fileName, matrixName, nx + ny, nx + nu);

    algorithm
      result.A := ABCD[1:nx, 1:nx];
      result.B := ABCD[1:nx, nx + 1:nx + nu];
      result.C := ABCD[nx + 1:nx + ny, 1:nx];
      result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
      Streams.print("State space record loaded from file: \"" +
        Modelica.Utilities.Files.fullPathName(fileName) + "\"");

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.Import.<b>fromFile</b>(fileName, matrixName)
</pre></blockquote>

<h4>Description</h4>
<p>
Reads and loads a state space system from a mat-file <tt>fileName</tt>. The file must contain the matrix [A, B; C, D] named matrixName and the integer nx representing the order of the system, i.e. the number of rows of the square matrix A.
</p>

<h4>Example</h4>
<blockquote><pre>

<b>algorithm</b>
  ss:=Modelica_LinearSystems2.StateSpace.Import.fromFile(\"stateSpace.mat\", \"ABCD\");
//  ss=StateSpace(
      A=[-1, 0, 0; 0, -2, 0; 0, 0, -3],
      B=[1; 1; 0],
      C=[1, 1, 1],
      D=[0])
</pre></blockquote>
</html>"));
    end fromFile;

    function fromModel
      "Generate a StateSpace data record by linearization of a model"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.StateSpace;
      import Simulator = DymolaCommands.SimulatorAPI;

      input String modelName "Name of the model"
       annotation(Dialog(__Dymola_translatedModel(translate=true)));
      input Real T_linearize = 0
        "Simulate until T_linearize and then linearize the model";
      input String fileName = "dslin" "Name of the result file";
      input String method = "Dassl" "Integration method";

    protected
      String fileName2 = fileName + ".mat"
        "Name of the result file with extension";
      Boolean OK1 = Simulator.simulateModel(
              problem=modelName,
              startTime=0,
              stopTime=T_linearize,
              method=method);
      Boolean OK2 = Simulator.importInitial("dsfinal.txt");
      Boolean OK3 = Simulator.linearizeModel(
              problem=modelName,
              resultFile=fileName,
              startTime=T_linearize,
              stopTime=T_linearize + 1,
              method=method);

    public
      output StateSpace result = StateSpace.Internal.read_dslin(fileName)
        "Model linearized at initial point";

    algorithm

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.Import.<b>fromModel</b>(modelName, T_linearize, fileName)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate a StateSpace data record by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation. The result of linearization is transformed into a StateSpace record.
</p>
<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
  Real T_linearize = 5;

<b>algorithm</b>
  ss = Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName, T_linearize);

// ss.A = [ 0.0,   1.0,    0.0,            0.0,      0.0,     0.0;
            0.0,   0.0,          -2.26,    0.08,     1.95,   -0.45;
            0.0,   0.0,           0.0,            1.0,      0.0,     0.0;
            0.0,   0.0,          -3.09,   -1.38,     7.70,   -3.01;
            0.0,   0.0,           0.0,            0.0,      0.0,     1.0;
            0.0,   0.0,          -6.47,    1.637,   -2.90,    1.29],

// ss.B=[0.0; 0.13; 0.0; -0.014; 0.0; -0.1],
// ss.C=identity(6),
// ss.D=[0; 0; 0; 0; 0; 0]
</pre></blockquote>
</html>", revisions="<html>
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
    end fromModel;

    annotation (Documentation(info="<html>
</html>"));
  end Import;

  encapsulated package Internal
    "Package of internal material of record StateSpace (for advanced users only)"
    extends Modelica.Icons.InternalPackage;

    import Modelica;
    import Modelica_LinearSystems2;

    encapsulated function assignOneOrTwoPoles
      "Algorithm to assign p (p = 1 or 2) eigenvalues"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Vectors;

      input Real F[:, size(F, 1)] "System matrix of order p=1 or p=2";
      input Real G[size(F, 1), :] "Control input matrix p rows";
      input Complex gamma[size(F, 1)];
      input Real tolerance=Modelica.Constants.eps;
      output Real K[:, size(F, 1)] "Feedback matrix p columns";

    protected
      Real Gamma[:, :];
      Integer rankGs;
      Real Fs[size(F, 1), size(F, 2)];
      Real Gs[size(G, 1), size(G, 2)];
      Real Gst[:, :]=transpose(G);
      Real Ks[:, size(F, 1)];
      Real c;
      Real s;
      Real r;
      Integer p=size(F, 1);
      Real sigmaG[:];

      Real V[size(G, 2), size(G, 2)];
      Real U[size(F, 1), size(F, 2)];

      Real u1[:];
      Real u2[:];
      Integer i;
      Complex system_ev[:];

    algorithm
      assert(size(F, 1) >= size(gamma, 1),
        "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F is of size ["
         + String(size(F, 1)) + "," + String(size(F, 1)) + "] and " + String(
        size(F, 1)) + " demanded assigned poles are expected. However, " +
        String(size(gamma, 1)) + " poles are given");
      //assert(not Modelica.Math.Matrices.isEqual(G,zeros(size(G,1),size(G,2)),tolerance),"A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix ");
      if size(gamma, 1) == 1 then
        assert(gamma[1].im == 0,
          "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F has size ["
           + String(size(F, 1)) + "," + String(size(F, 1)) +
          "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
           + String(gamma[1].im));
      elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
        assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im,
          "\nThe assigned pole pair given in function StateSpace.Internal.assignOneOrTwoPoles()"
          + "must be conjungated complex. However, the poles are"
          + "\npole1 = " + Complex.'String'(gamma[1])
          + "\npole2 = " + Complex.'String'(gamma[2])
          + ". \nTry\npole1 = " + Complex.'String'(gamma[1])
          + "\npole2 = " + Complex.'String'(
            Modelica.ComplexMath.conj(gamma[1])) + "\ninstead");
      end if;

      if not Modelica.Math.Matrices.isEqual(
              G,
              zeros(size(G, 1), size(G, 2)),
              tolerance) then
        if size(G, 2) == 1 then
          V := [1];
          if size(G, 1) == 1 then
            U := [1];
          else
            // Givens
            r := sqrt(G[1, 1]^2 + G[2, 1]^2);
            c := G[1, 1]/r;
            s := G[2, 1]/r;
            U := [c, s; -s, c];
          end if;
          Gs := U*G;

          rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
        else
          // size(G, 2)>1

          if size(G, 1) == 1 then
            // U=I, compute V by just one Householder transformation
            U := [1];
            u1 := cat(1, Vectors.householderVector(Gst[:, 1], cat(
                  1,
                  {1},
                  zeros(size(G, 2) - 1))));
            // Householder vector
            Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(
              Gst, u1);

            V := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1*
              u1);
            Gs := transpose(Gst);
            rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;

          else
            // systems with p==2 and m>1 are transformed by svd
            (sigmaG,U,V) := Modelica.Math.Matrices.singularValues(G);
            rankGs := 0;
            i := size(sigmaG, 1);
            while i > 0 loop
              if sigmaG[i] > 1e-10 then
                rankGs := i;
                i := 0;
              end if;
              i := i - 1;
            end while;
            Gs := zeros(p, size(G, 2));
            for i in 1:rankGs loop
              Gs[i, i] := sigmaG[i];
            end for;

          end if;
          V := transpose(V);
        end if;

        // check controllability
        assert(not Modelica.Math.Matrices.isEqual(
              Gs,
              zeros(size(Gs, 1), size(Gs, 2)),
              tolerance),
          "A subsystem in StateSpace.Internal.assignOneOrTwoPoles() is not controllable");

        Ks := fill(
              0,
              rankGs,
              size(F, 1));
        Fs := U*F*transpose(U);

        if size(F, 1) == 1 then
          Ks := matrix((Fs[1, 1] - gamma[1].re)/Gs[1, 1]);
        else
          if rankGs == size(F, 1) then

            // Gamma:= if size(F,1)==1 then [gamma[1].re] else [gamma[1].re, -(gamma[1].im)^2;1, gamma[2].re];
            //  Ks :=  Modelica_LinearSystems2.Math.Matrices.solve2(Gs, Fs - Gamma);
            //        Ks := [(Fs[1, 1] - gamma[1].re)/Gs[1, 1] - Gs[1, 2]*(Fs[2, 1] - 1)/Gs[1, 1]/Gs[2,2],
            //        (Fs[1, 2] + (gamma[1].im)^2)/Gs[1, 1] - Gs[1, 2]*(Fs[2, 2] - gamma[2].re)/Gs[1, 1]/Gs[2, 2];
            //        (Fs[2, 1] - 1)/Gs[2, 2],(Fs[2, 2] - gamma[2].re)/Gs[2, 2]];

            // since G1 is diagonal because of svd, Gs[1, 2] is zero

            Ks := [(Fs[1, 1] - gamma[1].re)/Gs[1, 1], (Fs[1, 2] + (gamma[1].im)
              ^2)/Gs[1, 1]; (Fs[2, 1] - 1)/Gs[2, 2], (Fs[2, 2] - gamma[2].re)/
              Gs[2, 2]];
          else

            Ks[1, 1] := (gamma[1].re + gamma[2].re - Fs[1, 1] - Fs[2, 2])/Gs[1,
              1];
            Ks[1, 2] := Ks[1, 1]*Fs[2, 2]/Fs[2, 1] + (Fs[1, 1]*Fs[2, 2] - Fs[1,
              2]*Fs[2, 1] - (gamma[1].re*gamma[2].re - gamma[1].im*gamma[2].im))
              /Fs[2, 1]/Gs[1, 1];
            Ks := -Ks;
          end if;
        end if;

        K := V[:, 1:size(Ks, 1)]*Ks*U;

      else
        if p == 1 then
          Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, " + String(F[1, 1]) + " remains and " + String(gamma[1].re) +
            " cannot be realized");
        else
          system_ev :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(F);
          Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, " + String(system_ev[1].re) + (if abs(system_ev[1].im) > 0 then
            " + " else " - ") + String(system_ev[1].im) + "j and " + String(
            system_ev[2].re) + (if abs(system_ev[2].im) > 0 then " + " else
            " - ") + String(system_ev[2].im) + "j remain and " + String(gamma[1].re)
             + (if abs(gamma[1].im) > 0 then (if gamma[1].im > 0 then " + "
             else " - " + String(gamma[1].im) + "j") else "" + " and ") +
            String(gamma[2].re) + (if abs(gamma[2].im) > 0 then (if gamma[2].im
             > 0 then " + " else " - " + String(gamma[2].im) + "j") else "") +
            " cannot be realized");
        end if;
        K := zeros(size(G, 2), size(F, 1));
      end if;

    end assignOneOrTwoPoles;

    encapsulated function assignOneOrTwoPoles_alpha
      "Algorithm to assign p (p = 1 or 2) eigenvalues"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Vectors;

      input Real F[:, size(F, 1)] "System matrix of order p=1 or p=2";
      input Real G[size(F, 1), :] "Control input matrix p rows";
      input Complex gamma[size(F, 1)];
      input Real tolerance=Modelica.Constants.eps;
      output Real K[:, size(F, 1)] "Feedback matrix p columns";

    protected
      Real Gamma[:, :];
      Integer rankGs;
      Real Fs[size(F, 1), size(F, 2)];
      Real Gs[size(G, 1), size(G, 2)];
      Real Gst[:, :]=transpose(G);
      Real Ks[:, size(F, 1)];
      Real c;
      Real s;
      Real r;

      Real V1[size(G, 2), size(G, 2)];
      Real V2[size(G, 2), size(G, 2)];
      Real V[size(G, 2), size(G, 2)];
      Real U[size(F, 1), size(F, 2)];

      Real u1[:];
      Real u2[:];

    algorithm
      assert(size(F, 1) >= size(gamma, 1),
        "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F is of size ["
         + String(size(F, 1)) + "," + String(size(F, 1)) + "] and " + String(
        size(F, 1)) + " demanded assigned poles are expected. However, " +
        String(size(gamma, 1)) + " poles are given");
      //assert(not Modelica.Math.Matrices.isEqual(G,zeros(size(G,1),size(G,2)),tolerance),"A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix ");
      if size(gamma, 1) == 1 then
        assert(gamma[1].im == 0,
          "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F has size ["
           + String(size(F, 1)) + "," + String(size(F, 1)) +
          "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
           + String(gamma[1].im));
      elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
        assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im,
          "\nThe assigned pole pair given in function StateSpace.Internal.assignOneOrTwoPoles()"
          + "must be conjungated complex. However, the poles are"
          + "\npole1 = " + Complex.'String'(gamma[1])
          + "\npole2 = " + Complex.'String'(gamma[2])
          + ". \nTry\npole1 = " + Complex.'String'(gamma[1])
          + "\npole2 = " + Complex.'String'(
            Modelica.ComplexMath.conj(gamma[1])) + "\ninstead");
      end if;

      if not Modelica.Math.Matrices.isEqual(
              G,
              zeros(size(G, 1), size(G, 2)),
              tolerance) then
        if size(G, 2) == 1 then
          V := [1];
          if size(G, 1) == 1 then
            U := [1];
          else
            // Givens
            r := sqrt(G[1, 1]^2 + G[2, 1]^2);
            c := G[1, 1]/r;
            s := G[2, 1]/r;
            U := [c, s; -s, c];
          end if;
          Gs := U*G;

          rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
        else
          // size(G, 2)>1

          if size(G, 1) == 1 then
            // U=I, compute V by just one Householder transformation
            U := [1];
            u1 := cat(1, Vectors.householderVector(Gst[:, 1], cat(
                  1,
                  {1},
                  zeros(size(G, 2) - 1))));
            // Householder vector
            Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(
              Gst, u1);

            V := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1*
              u1);
            Gs := transpose(Gst);
            rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;

          else
            //2xHH + Givens
            u1 := cat(1, Vectors.householderVector(Gst[:, 1], cat(
                  1,
                  {1},
                  zeros(size(G, 2) - 1))));
            // Householder vector1
            Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(
              Gst, u1);
            V1 := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1
              *u1);

            // if rank of G of a multi input system is equal to 1
            if Modelica.Math.Vectors.isEqual(
                    Gst[:, 2],
                    zeros(size(G, 2)),
                    tolerance) or Modelica.Math.Matrices.isEqual(
                    Gst[2:size(Gst, 1), :],
                    zeros(size(Gst, 1) - 1, size(Gst, 2)),
                    tolerance) then
              V := V1;
              rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
            else

              u2 := cat(
                    1,
                    zeros(1),
                    Vectors.householderVector(Gst[2:size(G, 2), 2], cat(
                      1,
                      {1},
                      zeros(size(G, 2) - 2))));
              // Householder vector2
              Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(
                Gst, u2);

              V2 := identity(size(G, 2)) - 2*matrix(u2)*transpose(matrix(u2))/(
                u1*u1);
              V := V2*V1;

            end if;

            Gs := transpose(Gst);

            rankGs := 0;
            for i in 1:2 loop
              if abs(Gs[i, i]) > tolerance then
                rankGs := rankGs + 1;
              end if;
            end for;
            // Givens rotation to transfotm Gs[1:2,1:2] to right upper triangle
            r := sqrt(Gs[1, 1]^2 + Gs[2, 1]^2);
            c := Gs[1, 1]/r;
            s := Gs[2, 1]/r;
            U := [c, s; -s, c];
            Gs := U*Gs;

          end if;
        end if;

        // check controllability
        assert(not Modelica.Math.Matrices.isEqual(
              Gs,
              zeros(size(Gs, 1), size(Gs, 2)),
              tolerance),
          "A subsystem in StateSpace.Internal.assignOneOrTwoPoles() is not controllable");

        Ks := fill(
              0,
              rankGs,
              size(F, 1));
        Fs := U*F*transpose(U);

        if size(F, 1) == 1 then
          Ks := matrix((Fs[1, 1] - gamma[1].re)/Gs[1, 1]);
        else
          if rankGs == size(F, 1) then

            // Gamma:= if size(F,1)==1 then [gamma[1].re] else [gamma[1].re, -(gamma[1].im)^2;1, gamma[2].re];
            //  Ks :=  Modelica_LinearSystems2.Math.Matrices.solve2(Gs, Fs - Gamma);
            Ks := [(Fs[1, 1] - gamma[1].re)/Gs[1, 1] - Gs[1, 2]*(Fs[2, 1] - 1)/
              Gs[1, 1]/Gs[2, 2], (Fs[1, 2] + (gamma[1].im)^2)/Gs[1, 1] - Gs[1,
              2]*(Fs[2, 2] - gamma[2].re)/Gs[1, 1]/Gs[2, 2]; (Fs[2, 1] - 1)/Gs[
              2, 2], (Fs[2, 2] - gamma[2].re)/Gs[2, 2]];
          else

            Ks[1, 1] := (gamma[1].re + gamma[2].re - Fs[1, 1] - Fs[2, 2])/Gs[1,
              1];
            Ks[1, 2] := Ks[1, 1]*Fs[2, 2]/Fs[2, 1] + (Fs[1, 1]*Fs[2, 2] - Fs[1,
              2]*Fs[2, 1] - (gamma[1].re*gamma[2].re - gamma[1].im*gamma[2].im))
              /Fs[2, 1]/Gs[1, 1];
            Ks := -Ks;
          end if;
        end if;

        //    K := transpose(V)*[Ks; zeros(size(G, 2) - rankGs, size(Ks, 2))]*(U);
        K := transpose(V[1:size(Ks, 1), :])*Ks*U;

      else
        Modelica.Utilities.Streams.print(
          "\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained");
        K := zeros(size(G, 2), size(F, 1));
      end if;

    end assignOneOrTwoPoles_alpha;

    function characterizeEigenvalue
      "Check stability, stabilizability, controllability, observability nad detectability of the single poles"

      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Eigenvalue;

      input StateSpace ss=Modelica_LinearSystems2.StateSpace(
              A=fill(0,0,0),
              B=fill(0,0,0),
              C=fill(0,0,0),
              D=fill(0,0,0)) "State space system";
      input Eigenvalue evin[:]
        "Eigenvalues or a pairs of conjugated complex pair";
      output Eigenvalue ev[size(ss.A, 1)];

    protected
      Real cPoles[:, 2] "Controllable poles";
      Real ncPoles[:, 2] "Uncontrollable poles";
      Real poles[size(ss.A, 1), 2] "Controllable and uncontrollable poles";
      Integer n_c;
      Integer nx=size(ss.A, 1);

      Eigenvalue cdummy[size(ss.A, 1)];
      Eigenvalue odummy[size(ss.A, 1)];
      StateSpace sst=StateSpace(
              A=transpose(ss.A),
              B=transpose(ss.C),
              C=transpose(ss.B),
              D=transpose(ss.D));
      Boolean equal;
      Integer ii;

      Real eps=1e-16*Modelica.Math.Matrices.norm(A=ss.A, p=1);
      Real factor_eps=1;
      Real absVector[nx];
      Integer indices[:];
      Integer indexMin;
      Real indexVector[:];
      Integer vv[:];

    algorithm
      for i in 1:nx loop
        ev[i].ev := evin[i].ev;
      end for;

      (cPoles,ncPoles,poles) := StateSpace.Internal.controllablePoles(ss);
      n_c := size(cPoles, 1);
      for i in 1:n_c loop
        cdummy[i].ev := cPoles[i, 1] + cPoles[i, 2]*j;
        cdummy[i].isControllable := true;
      end for;
      for i in 1:size(ncPoles, 1) loop
        cdummy[n_c + i].ev := ncPoles[i, 1] + ncPoles[i, 2]*j;
        cdummy[n_c + i].isControllable := false;
      end for;

      // controllable poles of the trnasposed system are the observable poles
      (cPoles,ncPoles,poles) := StateSpace.Internal.controllablePoles(sst);
      n_c := size(cPoles, 1);
      for i in 1:n_c loop
        odummy[i].ev := cPoles[i, 1] + cPoles[i, 2]*j;
        odummy[i].isObservable := true;
      end for;
      for i in 1:size(ncPoles, 1) loop
        odummy[n_c + i].ev := ncPoles[i, 1] + ncPoles[i, 2]*j;
        odummy[n_c + i].isObservable := false;
      end for;

      // using ev.imag as an flag to mark the eigenvalues
      for i in 1:nx loop
        absVector := fill(1e50, nx);
        for ii in 1:nx loop
          if not odummy[ii].imag then
            absVector[ii] :=Modelica.ComplexMath.'abs'(ev[i].ev - odummy[ii].ev);
          end if;
        end for;
        (absVector,indices) := Modelica.Math.Vectors.sort(absVector);
        indexMin := indices[1];
        ev[i].isObservable := odummy[indexMin].isObservable;
        odummy[indexMin].imag := true;
      end for;

      for i in 1:nx loop
        absVector := fill(1e50, nx);
        for ii in 1:nx loop
          if not cdummy[ii].imag then
            absVector[ii] :=Modelica.ComplexMath.'abs'(ev[i].ev - cdummy[ii].ev);
          end if;
        end for;
        (absVector,indices) := Modelica.Math.Vectors.sort(absVector);
        indexMin := indices[1];
        ev[i].isControllable := cdummy[indexMin].isControllable;
        cdummy[indexMin].imag := true;
      end for;

      for i in 1:nx loop
        ev[i] := Eigenvalue(
              ev[i].ev,
              ev[i].isControllable,
              ev[i].isObservable);
      end for;

    end characterizeEigenvalue;

    encapsulated function cntrHessenberg
      "Calculate the controllable part of a SISO system"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss "State space system";

      output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), 1],
        redeclare Real C[1, size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "controllable state space system";

    protected
      Integer nx=size(ss.A, 1);
      Real Ah1[nx, nx];
      Real bh1[nx];
      Real ch1[nx];
      Real u[:] "householder vector";
      Real Q[nx, nx];
      Real V[size(ss.A, 1), size(ss.A, 2)];
      Real tau[max(0,nx - 1)];
      Real Qc[:, :];
      Real svd[:];
      Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);
      Integer rankMinSys;
      Boolean isZero=false;

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      Ah1 := ss.A;
      bh1 := ss.B[:, 1];
      ch1 := ss.C[1, :];

      if Modelica.Math.Vectors.length(bh1) > 0 then

        if nx > 1 then

          // transform b->Qb = {*,0,...,0} and c->cQ, A->QAQ
          u := Vectors.householderVector(bh1, cat(
                1,
                {1},
                fill(0, nx - 1)));
          //householder vector to compute a housholder reflector S = I - 2*u*u'/u'*u
          Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);
          bh1 := Vectors.householderReflexion_e1(bh1, u);
          ch1 := Vectors.householderReflexion(ch1, u);

          (Ah1,V,tau) := Matrices.toUpperHessenberg(
                Ah1,
                1,
                nx);
          Q := Matrices.orthogonalQ(
                V,
                tau,
                1,
                nx);
          ch1 := ch1*Q;

        end if;

        rankMinSys := 1;
        while rankMinSys < nx and not isZero loop
          isZero := abs(Ah1[rankMinSys + 1, rankMinSys]) < normA*1e-10;
          rankMinSys := rankMinSys + 1;
        end while;

        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=Ah1,
              B=matrix(bh1),
              C=transpose(matrix(ch1)),
              D=ss.D,
              r=if isZero then rankMinSys - 1 else rankMinSys);

      end if;

      //equation

      //algorithm
    end cntrHessenberg;

    encapsulated function complexPoles
      "Generate a zeros-and-poles representation from state space representation"

      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Complex poles2[:]=fill(Complex(0, 0),
          StateSpace.Internal.numberOfPoles(ss));

    protected
      Integer nx=size(ss.A, 2);

      Complex zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";
      Complex zeros2[:];
      Complex poles[:] "eigenvalues of ss";
      Real eval[nx, 2];
      Real evec[nx, nx];

      Integer index[:]=fill(0, nx) "indices of zeros which are equal to poles";
      Integer i;
      Integer j;
      Integer k;
      Boolean h;
      Integer nzero;

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      h := false;
      zeros := StateSpace.Analysis.invariantZeros(ss);
      zeros2 := zeros;

      poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);

      for i in 1:size(zeros, 1) loop
        for j in 1:size(poles, 1) loop
          if zeros[i] == poles[j] then
            h := false;
            k := 1;
            while ((k < i) and (not h)) loop
              h := if (j == index[k]) then true else false;
              k := k + 1;
            end while;
            index[i] := if h then 0 else j;
          end if;
        end for;
      end for;

      j := 0;
      for i in 1:size(zeros, 1) loop
        if index[i] == 0 then
          j := j + 1;
          zeros2[j] := zeros[i];
        end if;
      end for;
      nzero := j;
      j := 0;
      for i in 1:size(poles, 1) loop
        h := false;
        k := 1;
        while (k <= size(zeros, 1) and (not h)) loop
          h := if i == index[k] then true else false;
          k := k + 1;
        end while;
        if not h then
          j := j + 1;
          poles2[j] := poles[i];

        end if;
      end for;

    end complexPoles;

    encapsulated function complexZeros
      "Calculate the zeros of the related transfer function"

      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Complex zeros2[:]=fill(Complex(0, 0),
          StateSpace.Internal.numberOfZeros(ss));

    protected
      Integer nx=size(ss.A, 2);

      Complex zeros[:];
      Complex poles[:];
      Real eval[nx, 2];
      Real evec[nx, nx];

      Integer index[:]=fill(0, nx) "Indices of zeros which are equal to poles";
      Integer i;
      Integer j;
      Integer k;
      Boolean h;
      Integer nzero;

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "Function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          "Function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      zeros := StateSpace.Analysis.invariantZeros(ss);

      poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);

      for i in 1:size(zeros, 1) loop
        for j in 1:size(poles, 1) loop
          if zeros[i] == poles[j] then
            h := false;
            k := 1;
            while ((k < i) and (not h)) loop
              h := if (j == index[k]) then true else false;
              k := k + 1;
            end while;
            index[i] := if h then 0 else j;
          end if;
        end for;
      end for;

      j := 0;
      for i in 1:size(zeros, 1) loop
        if index[i] == 0 then
          j := j + 1;
          zeros2[j] := zeros[i];
        end if;
      end for;

    end complexZeros;

    encapsulated function controllablePoles
      "Compute the controllable and uncontrollable poles of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal;

      input StateSpace ss=StateSpace(
              A=[-1],
              B=[1],
              C=[0],
              D=[0]) "State space system";

      output Real cPoles[:, 2] "controllable poles";
      output Real ncPoles[:, 2] "uncontrollable poles";
      output Real poles[size(ss.A, 1), 2]
        "controllable and uncontrollable poles";
    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssch(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "upper block controller Hessenberg form state space system";
      Boolean isControllable;

    algorithm
      if size(ss.B, 2) == 0 then
        poles := Modelica.Math.Matrices.eigenValues(ss.A);
        ncPoles := poles;
        cPoles := fill(
              0,
              0,
              2);
      else
        // build upper Hessenberg staircase to decomposite controllable/uncontrollable subspaces
        // The controllable part of A is in A[1:ssch.r, 1:ssch.r]
        (isControllable,ssch) := StateSpace.Internal.staircaseSVD(ss);
        if isControllable then
          poles := Modelica.Math.Matrices.eigenValues(ss.A);
          cPoles := poles;
          ncPoles := fill(
                0,
                0,
                2);
        else
          cPoles := Modelica.Math.Matrices.eigenValues(ssch.A[1:ssch.r, 1:ssch.r]);
          ncPoles := Modelica.Math.Matrices.eigenValues(ssch.A[ssch.r + 1:size(
            ss.A, 1), ssch.r + 1:size(ss.A, 1)]);
          poles := [cPoles; ncPoles];
        end if;
      end if;

      annotation (Documentation(info="<html>
The function uses the SVD based staircase algorithm to transform the state space representation into a similar state space
to separate the uncontrollable poles from the controllable poles.
</html>"));
    end controllablePoles;

    function damping "Frequencies and damping of state space system"
      extends Modelica.Icons.Function;

      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Complex eigenvalues[size(ss.A, 1)];
      output Real damp[size(ss.A, 1)];
      output Real frequency[size(ss.A, 1)];

    protected
      Integer n=size(ss.A, 1);
      Real pi=Modelica.Constants.pi;

    algorithm
      eigenvalues := StateSpace.Analysis.eigenValues(ss);
      for i in 1:n loop
        (frequency[i],damp[i]) :=Modelica_LinearSystems2.Math.ComplexAdvanced.frequency(eigenvalues[i]);
        frequency[i] := 2*pi*frequency[i];
      end for;

    end damping;

    encapsulated function dgreeOfRedSys
      "Calculate the controllable and observable part of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      output Integer degree_rs;

    protected
      StateSpace ssm1=
          Modelica_LinearSystems2.StateSpace.Internal.reducedCtrSystem(ss);

      StateSpace ss2=StateSpace(
              A=transpose(ssm1.A),
              B=transpose(ssm1.C),
              C=transpose(ssm1.B),
              D=ssm1.D);

      StateSpace ssm2=
          Modelica_LinearSystems2.StateSpace.Internal.reducedCtrSystem(ss2);

    algorithm
      degree_rs := size(ssm2.A, 1);

    end dgreeOfRedSys;

    encapsulated function householder "Householder transformation"
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss;
      input Real v[size(ss.A, 1)];

      output StateSpace ssh;

    protected
      Real Ah[size(ss.A, 1), size(ss.A, 2)];
      Real Bh[size(ss.B, 1), size(ss.B, 2)];
      Real Ch[size(ss.C, 1), size(ss.C, 2)];
    algorithm
      Ah := Matrices.householderSimilarityTransformation(ss.A, v);
      Bh := Matrices.householderReflexion(ss.B, v);
      Ch := transpose(Matrices.householderReflexion(transpose(ss.C), v));

      ssh := StateSpace(
            A=Ah,
            B=Bh,
            C=Ch,
            D=ss.D);

    end householder;

    encapsulated function invariantZerosWithRealMatrix
      "Compute invariant zeros of linear state space system (system given by A,B,C,D matrices)"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Complex;

      input Real A[:,size(A,1)] "A-matrix of linear state space system";
      input Real B[size(A,1),:] "B-matrix of linear state space system";
      input Real C[:,size(A,1)] "C-matrix of linear state space system";
      input Real D[size(C,1),size(B,2)] "D-matrix of linear state space system";

      output Real InvariantZeros[:,2]
        "Finite, invariant zeros of linear state space system; size(Zeros,1) <= size(A,1); Zeros[:,1]: Real part, Zeros[:,2]: Imaginary part";

    protected
      Integer n;
      Integer m;
      Integer p;
      Real Ar[:, :];
      Real Br[:, :];
      Real Cr[:, :];
      Real Dr[:, :];

      Real Af[:, :];
      Real Bf[:, :];
      Real AfBf[:, :];

      Real V2[:, :];
      Real Vf[:, :];
      Real R[:, :];

      Integer na;
      Real alphaReal[:];
      Real alphaImag[:];
      Real beta[:];
      Integer info;
      Real zerosMax;
      Real absZero[:];

      Integer j;
    algorithm
      if min(size(B)) == 0 or min(size(C)) == 0 then
        // no inputs or no outputs
        InvariantZeros := fill(0.0,0,2);
      else
        // Reduce system
        (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(A,B,C,D);
        if n > 0 then
          (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(
                transpose(Ar),
                transpose(Cr),
                transpose(Br),
                transpose(Dr));
        end if;
        if n == 0 then
           InvariantZeros := fill(0.0,0,2);
        else
          (,R,,V2) := Matrices.QR(Matrices.fliplr(transpose([Cr, Dr])));
          Vf := Matrices.fliplr(V2);
          AfBf := [Ar, Br]*Vf;
          Af := AfBf[:, 1:size(Ar, 2)];
          Bf := Vf[1:size(Ar, 1), 1:size(Ar, 2)];

          (alphaReal,alphaImag,beta,info) :=
             Modelica_LinearSystems2.Math.Matrices.LAPACK.dggev_eigenValues(Af,Bf);
          assert(info == 0,
            "Failed to compute invariant zeros with function invariantZerosWithRealMatrix(..)");

          InvariantZeros := fill(0.0, size(beta, 1),2);

          /* The pencil (Af,Bf) has n eigenvalues, since the transformation to this
             form is done so that Bf is regular. Therefore, the generalized eigenvalues
             represented by alpha[i]/beta[i] have the property that beta[i] cannot be zero
             and a division by beta[i] is uncritical.

             |alpha| / beta <= zerosMax
             if |alpha| <= beta*zerosMax then
                zero is used
             else
                assumed that zero is at infinity (i.e. it is ignored)
          */
          j := 0;
          zerosMax := 1.0e4*Modelica.Math.Matrices.norm([Af, Bf], p=1);
          absZero := Modelica_LinearSystems2.StateSpace.Internal.absComplexVector(alphaReal, alphaImag);
          for i in 1:size(beta, 1) loop
             if absZero[i] <= beta[i]*zerosMax then
                j := j + 1;
                InvariantZeros[j,1] := alphaReal[i]/beta[i];
                InvariantZeros[j,2] := alphaImag[i]/beta[i];
             end if;
          end for;

          if j == 0 then
             InvariantZeros := fill(0.0, 0,2);
          else
             InvariantZeros := InvariantZeros[1:j,:];
          end if;
        end if;
      end if;
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
InvariantZeros = StateSpace.Internal.<b>invariantZerosWithRealMatrix</b>(A,B,C,D)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
     <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre></blockquote>
<p>
The invariant zeros of this system are defined as the variables
s  that make the Rosenbrock matrix of the system
</p>
<pre>
    | s<b>I-A</b>   <b>-B</b> |
    |           |
    | <b>C</b>       <b>D</b> |
</pre>
<p>
singular.
</p>
<p>
This function applies the algorithm described in [1] where the system (<b>A</b>, <b>B</b>, <b>C</b>, <b>D</b>) is reduced to a new system (<b>A</b>r, <b>B</b>r <b>C</b>r, <b>D</b>r) with the same zeros and with <b>D</b>r of full rank.
</p>

<p>
The zeros are returned as a matrix InvariantZeros[:,2] where InvariantZeros[i,1] is the real and InvariantZeros[i,2] is the imaginary part
of the complex zero i.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[1, 1, 1;0, 1, 1;0,0,1],
    B=[1;0;1],
    C=[0,1,1],
    D=[0]);

  Complex zeros[:];

<b>algorithm</b>
  zeros := Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);
// zeros = {1, 0}
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Emami-Naeini, A. and Van Dooren, P. (1982):</dt>
<dd> <b>Computation of Zeros of Linear Multivariable Systems</b>.
     Automatica, 18, pp. 415-430.<br>&nbsp;</dd>
</dl>
</html>", revisions="<html>
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
    end invariantZerosWithRealMatrix;

    encapsulated function invariantZeros2
      "Compute invariant zeros of linear SISO state space system with a generalized system matrix [A, B, C, D] which is of upper Hessenberg form"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "Linear system in state space form";
      output Complex Zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
      Integer nx=size(ss.A, 1) "Number of states";
      Integer nu=size(ss.B, 2) "Number of inputs";
      Integer ny=size(ss.C, 1) "Number of outputs";
      Integer na=nx + nu;
      Real A[nx + ny, nx + nu]=[ss.A, ss.B; ss.C, ss.D];
      Real B[nx + ny, nx + nu]=[identity(nx), zeros(nx, nu); zeros(ny, nx + nu)];
      Real alphaReal[na];
      Real alphaImag[na];
      Real beta[na];
      Integer info;
      Real beta_small=100*Modelica.Constants.eps;
      Integer nZeros;
      Complex z[size(ss.A, 1)];
      Integer j;
      Real normB=max(beta_small,Modelica.Math.Matrices.norm(ss.B, p=1));
    algorithm
      assert(nu == ny, "Function invariantZeros requires currently that the number of
inputs (= " + String(nu) + ") = number of outputs (= " + String(ny) + ")
This condition is however not fulfilled");

      // Compute zeros

      (alphaReal,alphaImag,beta,info) :=
        Matrices.generalizedEigenvaluesTriangular(A, B);

      assert(info == 0,
        "Failed to compute invariant zeros with function invariantZeros(..)");

      // If beta[i] is zero, then zero i is infinite.
      j := 1;
      for i in 1:na loop

        if beta[i] >= normB*1e-10 then
          // finite eigenvalue
          z[j].re := if abs(alphaReal[i]) >= normB*1e-12 then alphaReal[i]/beta[
            i] else 0;
          z[j].im := if abs(alphaImag[i]) >= normB*1e-12 then alphaImag[i]/beta[
            i] else 0;
          j := j + 1;
        end if;
      end for;
      nZeros := j - 1;
      Zeros := z[1:nZeros];
      annotation (Documentation(info="<html>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
z that make the following matrix singular:
</p>
<pre>
    | <b>A</b> <b>B</b> |     | <b>I</b> <b>0</b> |
    |     | - z*|     |
    | <b>C</b> <b>D</b> |     | <b>0</b> <b>0</b> |
</pre>
<p>
where <b>I</b> is the identity matrix of the same size as <b>A</b>
and <b>0</b> are zero matrices of appropriate dimensions.
</p>
<p>
Unlike to function StateSpace.Analysis.invariantZeros for general systems, it is
assumned in StateSpace.Analysis.invariantZeros that the generalized system matrix
[<b>A</b>, <b>B</b>; <b>C</b>, <b>D</b>] has upper Hessenberg form. Especially for SISO system this is
achieved when <b>A</b> is of upper Hessenberg form and [1, n] matrix <b>C</b> is of form
<b>C</b> = k*[0, 0, ..., 0, 1].
<p>
The function uses the LAPACK routine DHGEQZ. Look at <b>Modelica_LinearSystems2.Math.Matrices.LAPACK.dhgeqz</b> for details.
<p>
The advantage of this function in comparison to the general invariantZeros function
is the lower computatioal effort bacause systems with arbitrary system functions are first transformed
into an upper Hessenberg form system.
<p>
This function is used in fromStateSpace transformation functions which use Hessenberg form systems anyway.
</p>
<p>
Currently, there is the restriction that the number of
inputs and the number of outputs must be identical. Other systems
have to be treated like p*q SISO systems where p is the number of putputs and q the number of inputs of the MIMO system.
</p>
</html>"));
    end invariantZeros2;

    encapsulated function invariantZerosHessenberg
      "Fast version to calculate the system zeros of a SISO system with D=[0] and A has upper Hessenberg form, delivered by StateSpace.reduceSystem"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "Linear system in state space form";
      output Complex Zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
      Integer nx=size(ss.A, 1) "Number of states";
      Integer nu=size(ss.B, 2) "Number of inputs";
      Integer ny=size(ss.C, 1) "Number of outputs";
      Real Ah[nx, nx]=ss.A;
      Real eps=100*Modelica.Constants.eps;
      Integer k;
      Boolean h;

    algorithm
      assert(nu == 1, "Function invariantZeros2 requires currently a SISO-input system.\n
This condition is however not fulfilled because the number of inputs is nu = "
         + String(nu));
      assert(ny == 1, "Function invariantZeros2 requires currently a SISO-input system.\n
This condition is however not fulfilled because the number of outputs is ny = "
         + String(ny));

      h := true;
      k := nx + 1;

      if size(ss.B, 2) > 0 then
        while k >= 1 and h loop
          k := k - 1;
          if abs(ss.B[k, 1]) >= eps then

            h := false;
          end if;

        end while;

        Zeros := fill(Complex(0), k - 1);

        if k > 1 then
          Ah[:, k - 1] := ss.A[:, k - 1] - ss.A[k, k - 1]/ss.B[k, 1]*ss.B[:, 1];

          //    Zeros := Complex.eigenValues(Ah[1:k - 1, 1:k - 1]);
          Zeros :=Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.eigenValues_dhseqr(Ah[1:k - 1, 1:k - 1]);

          for i in 1:k - 1 loop
            if Modelica.ComplexMath.'abs'(Zeros[i]) <
              Modelica.Math.Matrices.norm(Ah[1:k - 1, 1:k - 1], p=1)*1e-12 then
              Zeros[i] := Complex(0);
            end if;

          end for;

        end if;
      else
        Zeros := fill(Complex(0, 0), 0);
      end if;

      annotation (Documentation(info="<html>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
z that make the following matrix singular:
</p>
<pre>
    | <b>A</b> <b>B</b> |     | <b>I</b> <b>0</b> |
    |     | - z*|     |
    | <b>C</b> <b>D</b> |     | <b>0</b> <b>0</b> |
</pre>
<p>
where <b>I</b> is the identity matrix of the same size as <b>A</b>
and <b>0</b> are zero matrices of appropriate dimensions.
</p>
<p>
Currently, there is the restriction that the number of
inputs and the number of outputs must be identical.
</p>
</html>"));
    end invariantZerosHessenberg;

    encapsulated function isSISO
      "To check a state space system to be SISO (or not)"

      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean isSISO "True, if state space system is SISO";
    algorithm
      isSISO := size(ss.B, 2) == 1 and size(ss.C, 1) == 1;
      annotation (Documentation(info="<html>
</html>"));
    end isSISO;

    encapsulated function isControllableAndObservableSISO
      "To check whether a SISO system is controllable and observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssm1=
          StateSpace.Internal.reducedCtrSystem(ss);
      Integer nx=size(ss.A, 1);
      Integer rankQ=ssm1.r;
      StateSpace ss2=StateSpace(
              A=transpose(ssm1.A[nx - rankQ + 1:nx, nx - rankQ + 1:nx]),
              B=transpose(ssm1.C[:, nx - rankQ + 1:nx]),
              C=transpose(ssm1.B[nx - rankQ + 1:nx, :]),
              D=ssm1.D);
      Integer nx2=ssm1.r;
      Modelica_LinearSystems2.Internal.StateSpaceR ssm2=
          StateSpace.Internal.reducedCtrSystem(ss2);
    public
      output Boolean controllableAndObservable;
    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;
      controllableAndObservable := size(ss.A, 1) == ssm2.r;

    end isControllableAndObservableSISO;

    encapsulated function isControllableSISO
      "To check a SISO system whether it is controllable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssm=
          StateSpace.Internal.reducedCtrSystem(ss);

    public
      output Boolean controllable;
    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;
      controllable := size(ss.A, 1) == ssm.r;

      annotation (Documentation(info="<html>
This function is to calculate whether a SISO state space system is controllable or not. Therefore,
it is transformed to lower controller Hessenberg form
<blockquote><pre>
               | *    *     0   ...  0 |               | 0 |
               | .    .     .    .   . |               | . |
 <b>Q</b>*<b>A</b>*<b>Q</b> ' = <b>H</b> = | *   ...   ...   *   0 |,    <b>Q</b>*<b>b</b> = <b>q</b> = | . |,   <b>c</b>*<b>Q</b> = ( *, ..., * )
               | *   ...   ...   *   * |               | 0 |
               | *   ...   ...   *   * |               | * |

</pre>
</blockquote>
Note, that
<blockquote><pre>
                   n-1                        n-1
rank(<b>b</b>, <b>A</b>*<b>b</b>, ..., <b>A</b>  *<b>b</b>) = rank(<b>q</b>, <b>H</b>*<b>q</b>, ..., <b>H  </b>*<b>q</b>)
</pre>
</blockquote>
and that
<blockquote><pre>
                 n-1
 (<b>q</b>, <b>H</b>*<b>q</b>, ..., <b>H</b>  *<b>q</b>)
</pre>
</blockquote>
is a lower triangular matrix and has full rank if and only if none of the elements in the diagonal is zero. That is, that neither qn or hi,i+1,   i = 1,..., n-1   may be zero.
</html>"));
    end isControllableSISO;

    encapsulated function isControllableMIMO
      "To check a MIMO system whether it is controllable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

      output Boolean controllable;
    algorithm
      assert(method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD or method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.QR, "\nMethods for staircase algorithm are QR factorization or singular value decomposition. Therefore,
the variable \"method\" in \"Modelica_LinearSystems2.StateSpace.Internal.isControllableMIMO\" has to be qr or svd but is method = " + String(method));
      if min(size(ss.B)) == 0 then
        controllable := false;
      else
        if method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.QR then
          controllable := StateSpace.Internal.staircaseQR(ss);
        else
          controllable := StateSpace.Internal.staircaseSVD(ss);
        end if;
      end if;
      annotation (Documentation(info="<html>
</html>"));
    end isControllableMIMO;

    encapsulated function isDetectableSISO
      "To check whether a SISO system is detectable"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Boolean detectable;
    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssm=
          StateSpace.Internal.cntrHessenberg(
          StateSpace.Internal.transposeStateSpace(ss));
      Complex evd[:]=fill(Complex(0), size(ss.A, 1) - ssm.r);
    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      detectable := true;
      if size(ss.A, 1) > ssm.r then
        evd :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ssm.A[ssm.r + 1:size(ss.A, 1), ssm.r + 1:size(ss.A, 1)]);
        for i1 in 1:size(evd, 1) loop
          detectable := detectable and evd[i1].re < 0;
        end for;
      end if;

      annotation (Documentation(info="<html>
<p>
This function checks whether a SISO state space system is detectable or not.
</p>
<p>
A system is detectable for the continuous-time case if all of the unobservable eigenvalues have neagtive real part
or for the discrete-time case if all of the unobservable eigenvalues are in the complex unit circle respectively.
Hence, a oberservable system is always detectable of course.
</p>
<p>
As observability is a dual concept of controllability, the concept of detectability is dual to stabilizability, that is,
a system is detectable if the pair (<b>A</b>', <b>C</b>') is stabilizable. Therefore, the same algorithm to check stabilizability
are applied to the dual pair (<b>A</b>', <b>C</b>') of the system:
</p>
<p>
To check stabilizability (see Modelica_LinearSystems2.StateSpace.Analysis.isStabilizable) , ths system is transformed to to upper controller Hessenberg form
</p>
<blockquote><pre>
              | *   *   ...   ...    * |               | * |
              | *   *   ...   ...    * |               | 0 |
<b>Q</b>*<b>A</b>*<b>Q</b> ' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>b</b> = <b>q</b> = | . |,   <b>c</b>*<b>Q</b> = ( *, ..., * )
              | .   .    .     .     . |               | . |
              | 0  ...   0     *     * |               | 0 |
</pre>
</blockquote>
<p>
The system can be partitioned to
</p>
<blockquote><pre>
<b>H</b>=[<b>H</b>11,<b>H</b>12; <b>H</b>21, <b>H</b>22], <b>q</b>=[<b>q</b>1;<b>0</b>],
</pre></blockquote>
<p>
where the pair (<b>H</b>11, <b>q</b>1) contains the controllable part of the system, that is, rank(<b>H</b>) = rank(<b>H</b>11). For
stabilizability the <b>H</b>22 has to be stable.
</p>
</html>"));
    end isDetectableSISO;

    encapsulated function isDetectableMIMO
      "Check whether a MIMO system is detectable"

      import Modelica;
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean detectable;

    protected
      StateSpace sst=StateSpace.Internal.transposeStateSpace(ss);

      Complex evnd[:] "Complex vector of uncontrollable poles";
      Real dPoles[:, 2] "Controllable poles";
      Real ndPoles[:, 2] "Uncontrollable poles";
      Real poles[size(ss.A, 1), 2] "Controllable and uncontrollable poles";

    algorithm
      (dPoles,ndPoles,poles) := StateSpace.Internal.controllablePoles(sst);
      evnd := fill(Complex(0), size(ndPoles, 1));
      for i1 in 1:size(ndPoles, 1) loop
        evnd[i1] := ndPoles[i1, 1] + j*ndPoles[i1, 2];
      end for;

      detectable := true;

      if size(sst.A, 1) == size(dPoles, 1) then
        detectable := true;
      else
        for i1 in 1:size(ndPoles, 1) loop
          detectable := detectable and ndPoles[i1, 1] < 0;
        end for;
      end if;

      annotation (Documentation(info="<html>
This function checks whether a MIMO state space system is detectable or not.
<p>
A system is detectable for the continuous-time case if all of the unobservable eigenvalues have negative real part
or for the discrete-time case if all of the unobservable eigenvalues are in the complex unit circle respectively.
Hence, a observable system is always detectable of course.
<p>
To check detectability, staircase algorithm is used to separate the observable subspace from the unobservable subspace.
The unobservable poles are checked to be stable.
</html>"));
    end isDetectableMIMO;

    encapsulated function isObservableSISO
      "To check whether a SISO system is observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

    protected
      StateSpace ss2=StateSpace.Internal.transposeStateSpace(ss);

      Modelica_LinearSystems2.Internal.StateSpaceR ssm2=
          StateSpace.Internal.reducedCtrSystem(ss2);
    public
      output Boolean observable;
    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;
      observable := size(ss.A, 1) == ssm2.r;

      annotation (Documentation(info="<html>
<p>
This function is to calculate whether a SISO state space system is observable or not. Therefore, the dual System (A', c', b', d')
it is transformed to upper observer Hessenberg form
</p>
<blockquote><pre>
              | *   *   ...   ...    * |             | * |
              | *   *   ...   ...    * |             | . |
<b>Q</b>*<b>A'</b>*<b>Q</b>' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>c'</b> =  | . |,   <b>b'</b>*<b>Q</b> = <b>q</b> = ( 0, ..., 0, * )
              | .   .    .     .     . |             | * |
              | 0  ...   0     *     * |             | * |
</pre></blockquote>
<p>
Note, that
</p>
<blockquote><pre>
rank(<b>c'</b>; <b>c'</b>*<b>A'</b>; ...; <b>c'</b>*<b>A'</b><sup><big>(n-1)</big></sup>) = rank(<b>q</b>; <b>q</b>*<b>H</b>; ...; <b>q</b>*<b>H</b><sup><big>(n-1)</big></sup>)
</pre></blockquote>
<p>
and that
</p>
<blockquote><pre>
(<b>q</b>; <b>H</b>*<b>q</b>; ...; <b>q</b>*<b>H</b><sup><big>(n-1)</big></sup>)
</pre></blockquote>
<p>
is a lower triangular matrix and has full rank if and only if none of the elements in
the diagonal is zero. That is, that neither qn or hi,i-1,   i = 2,...,&nbsp;n   may be zero.
</p>
</html>"));
    end isObservableSISO;

    encapsulated function isObservableMIMO
      "To check a MIMO system whether it is observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      input Modelica_LinearSystems2.Utilities.Types.StaircaseMethod method=Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD;

    protected
      StateSpace ss2=StateSpace.Internal.transposeStateSpace(ss);

    public
      output Boolean observable;
    algorithm
      assert(method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.SVD or method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.QR, "\nMethods for staircase algorithm are QR factorization or singular value decomposition. Therefore,
the variable \"method\" in \"Modelica_LinearSystems2.StateSpace.Internal.isControllableMIMO\" has to be qr or svd but is method = " + String(method));
      if min(size(ss.C)) == 0 then
        observable := false;
      else
        if method == Modelica_LinearSystems2.Utilities.Types.StaircaseMethod.QR then
          observable := StateSpace.Internal.staircaseQR(ss2);
        else
          observable := StateSpace.Internal.staircaseSVD(ss2);
        end if;
      end if;

      annotation (Documentation(info="<html>
</html>"));
    end isObservableMIMO;

    encapsulated function isStabilizableSISO
      "To check whether a SISO system is stabliziable"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";
      output Boolean stabilizable;
    protected
      Modelica_LinearSystems2.Internal.StateSpaceR ssm=
        StateSpace.Internal.cntrHessenberg(ss);
      Complex evd[:]=fill(Complex(0), size(ss.A, 1) - ssm.r);

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;
      stabilizable := true;

      if size(ss.A, 1) > ssm.r then
        evd :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ssm.A[ssm.r + 1:size(ss.A, 1), ssm.r + 1:size(ss.A, 1)]);
        for i1 in 1:size(evd, 1) loop
          stabilizable := stabilizable and evd[i1].re < 0;
        end for;
      end if;

      annotation (Documentation(info="<html>
<p>
This function checks whether a SISO state space system is stabilizable or not.
</p>
<p>
A system is stabilizable for the continuous-time case if all of the uncontrollable eigenvalues have neagtive real part
or for the discrete-time case if all of the uncontrollable eigenvalues are in the complex unit circle respectively.
Hence, a controllable system is always stabilizable of course.
</p>
<p>
To check stabilizability, ths system is transformed to to upper controller Hessenberg form
</p>
<blockquote><pre>
              | *   *   ...   ...    * |               | * |
              | *   *   ...   ...    * |               | 0 |
<b>Q</b>*<b>A</b>*<b>Q</b> ' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>b</b> = <b>q</b> = | . |,   <b>c</b>*<b>Q</b> = ( *, ..., * )
              | .   .    .     .     . |               | . |
              | 0  ...   0     *     * |               | 0 |
</pre></blockquote>
<p>
The system can be partitioned to
</p>
<blockquote><pre>
<b>H</b>=[<b>H</b>11,<b>H</b>12; <b>H</b>21, <b>H</b>22], <b>q</b>=[<b>q</b>1;<b>0</b>],
</pre></blockquote>
<p>
where the pair (<b>H</b>11, <b>q</b>1) contains the controllable part of the system, that is, rank(<b>H</b>) = rank(<b>H</b>11). For
stabilizability the <b>H</b>22 has to be stable.
</p>
</html>"));
    end isStabilizableSISO;

    encapsulated function isStabilizableMIMO
      "To check whether a MIMO system is stabliziable"

      import Modelica;
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean stabilizable;

    protected
      Complex evnc[:] "complex vector of uncontrollable poles";
      Real cPoles[:, 2] "controllable poles";
      Real ncPoles[:, 2] "uncontrollable poles";
      Real poles[size(ss.A, 1), 2] "controllable and uncontrollable poles";

    algorithm
      (cPoles,ncPoles,poles) := StateSpace.Internal.controllablePoles(ss);
      evnc := fill(Complex(0), size(ncPoles, 1));
      for i1 in 1:size(ncPoles, 1) loop
        evnc[i1] := ncPoles[i1, 1] + j*ncPoles[i1, 2];
      end for;

      stabilizable := true;

      if size(ss.A, 1) > size(cPoles, 1) then
        for i1 in 1:size(ncPoles, 1) loop
          stabilizable := stabilizable and ncPoles[i1, 1] < 0;
        end for;
      end if;

      annotation (Documentation(info="<html>
This function checks whether a MIMO state space system is stabilizable or not.
<p>
A system is stabilizable for the continuous-time case if all of the uncontrollable eigenvalues have neagtive real part
or for the discrete-time case if all of the uncontrollable eigenvalues are in the complex unit circle respectively.
Hence, a controllable system is always stabilizable of course.
<p>
To check stabilizability, staircase algorithm is used to separate the controllable subspace from the uncontrollable subspace.
The uncontrollable poles are checked to to stable.
</html>"));
    end isStabilizableMIMO;

    encapsulated function numberOfPoles
      "Calculate the number of poles of the related transfer function"

      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Integer numberOfPoles=StateSpace.Internal.numberOfPolesAndZeros(ss);

    algorithm
    end numberOfPoles;

    encapsulated function numberOfPolesAndZeros
      "Calculate the number poles and of zeros of the related transfer function"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

      input StateSpace ss "State space system";
      output Integer numberOfPoles;
      output Integer numberOfZeros;

    protected
      Integer nx=size(ss.A, 2);
      Complex zeros[:];
      Complex zeros2[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";
      Complex poles[:];
      Complex poles2[:] "eigenvalues of ss";
      Real eval[nx, 2];
      Real evec[nx, nx];

      Integer index[:]=fill(0, nx) "indices of zeros which are equal to poles";
      Integer i;
      Integer j;
      Integer k;
      Boolean h;
      Integer nzero;

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      zeros := StateSpace.Analysis.invariantZeros(ss);
      zeros2 := zeros;

      poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);
      poles2 := poles;

      //Reduce terms which are in nominator as well as in denominator
      for i in 1:size(zeros, 1) loop
        for j in 1:size(poles, 1) loop
          if zeros[i] == poles[j] then
            h := false;
            k := 1;
            while ((k < i) and (not h)) loop
              h := if (j == index[k]) then true else false;
              k := k + 1;
            end while;
            index[i] := if h then 0 else j;
          end if;
        end for;
      end for;

      j := 0;
      for i in 1:size(zeros, 1) loop
        if index[i] == 0 then
          j := j + 1;
          zeros2[j] := zeros[i];
        end if;
      end for;
      nzero := j;
      j := 0;
      for i in 1:size(poles, 1) loop
        h := false;
        k := 1;
        while (k <= size(zeros, 1) and (not h)) loop
          h := if i == index[k] then true else false;
          k := k + 1;
        end while;
        if not h then
          j := j + 1;
          poles2[j] := poles[i];

        end if;
      end for;

      numberOfPoles := nx - size(zeros, 1) + nzero;
      numberOfZeros := nzero;

    end numberOfPolesAndZeros;

    encapsulated function numberOfZeros
      "Calculate the number of zeros of the related transfer function"

      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "StateSpace object";

      output Integer numberOfZeros;

    protected
      Integer nop;
      Integer noz;

    algorithm
      (nop,noz) := StateSpace.Internal.numberOfPolesAndZeros(ss);
      numberOfZeros := noz;
    end numberOfZeros;

    encapsulated function partialGain "Algorithm for partial gain"
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input Real H[:, size(H, 1)] "Upper Hessenberg matrix";
      input Real b[size(H, 1)];
      output Real result;
    protected
      Real Hh[:, :]=H;
      Real bh[:]=b;
      Integer q=size(H, 1);
    algorithm

      (Hh,bh) := StateSpace.Internal.trianUpperHess(Hh, bh);
      result := bh[q]/Hh[q, q];

    end partialGain;

    encapsulated function polesAndZeros
      "Generate poles and zeros from state space representation"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.PolesAndZeros;
      import Modelica_LinearSystems2.Internal;

      input StateSpace ss "State space system";
      input StateSpace ssm=
          Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(
          ss);
      output Internal.PolesAndZeros pz(
        redeclare Real p_real[size(ssm.A, 1)],
        redeclare Real p_im[size(ssm.A, 1)],
        redeclare Real z_real[size(StateSpace.Analysis.invariantZeros(ssm), 1)],
        redeclare Real z_im[size(StateSpace.Analysis.invariantZeros(ssm), 1)]);

    protected
      Complex poles[:]=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ssm.A);
      Complex zeros[:]=StateSpace.Analysis.invariantZeros(ssm);

    algorithm
      pz.p_real := poles[:].re;
      pz.p_im := poles[:].im;
      pz.z_real := zeros[:].re;
      pz.z_im := zeros[:].im;
      pz.norz_p := Internal.numberOfRealZeros(poles);
      pz.norz_z := Internal.numberOfRealZeros(zeros);

    end polesAndZeros;

    encapsulated function readLength_nu
      "Read the number of inputs nu of a state space system from a file"

      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

      input String fileName = "ss_siso.mat"
        "Name of the state space system data file" annotation (Dialog(
            loadSelector(filter="MAT files (*.mat);; All files (*.*)", caption=
                "State space system data file")));
      input String matrixName = "ABCD" "Name of the state space system matrix"
        annotation (Dialog);

      output Integer nu;
    protected
      Integer xuy[3] = ReadSystemDimension(fileName, matrixName);
    algorithm
      nu := xuy[2];
    end readLength_nu;

    encapsulated function readLength_nx
      "Read the order nx of a state space system from a file"

      import Modelica.Utilities.Streams;

      input String fileName = "ss_siso.mat"
        "Name of the state space system data file" annotation (Dialog(
            loadSelector(filter="MAT files (*.mat);; All files (*.*)", caption=
                "state space system data file")));
      output Integer nx;

    algorithm
      nx :=integer(scalar(Streams.readRealMatrix(fileName, "nx", 1, 1)));
    end readLength_nx;

    encapsulated function readLength_ny
      "Read the number of outputs ny of a state space system from a file"

      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

      input String fileName = "ss_siso.mat"
        "Name of the state space system data file" annotation (Dialog(
            loadSelector(filter="MAT files (*.mat);; All files (*.*)", caption=
                "state space system data file")));
      input String matrixName = "ABCD" "Name of the state space system matrix"
        annotation (Dialog);

      output Integer ny;
    protected
      Integer xuy[3] = ReadSystemDimension(fileName, matrixName);
    algorithm
      ny := xuy[3];
    end readLength_ny;

    encapsulated function reducedCtrSystem2
      "Calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss "State space system";
      input Real eps=0;

    protected
      StateSpace sst(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "Transformed state space system";

      Integer nx=size(ss.A, 1);
      Real Ah1[size(ss.A, 1), size(ss.A, 2)];
      Real bh1[size(ss.A, 1)];
      Real ch1[size(ss.A, 1)];

      Real u[:] "householder vector";
      Real cpoles[:, 2]=
          Modelica_LinearSystems2.StateSpace.Internal.controllablePoles(ss);
      Integer rankQc=size(cpoles, 1);

      Integer rankQc2;
      Real Qc2[nx, nx];
      Real sigma[:];
      Real eps2;

      Real Ah2[rankQc, rankQc];
      Real bh2[rankQc];
      Real ch2[rankQc];

      Integer ll;
      Integer r;

      Boolean h;

    public
      output StateSpace ssm1(
        redeclare Real A[rankQc, rankQc],
        redeclare Real B[rankQc, 1],
        redeclare Real C[1, rankQc],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "controllable state space system";

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      Ah1 := ss.A;
      bh1 := ss.B[:, 1];
      ch1 := ss.C[1, :];

      if nx > 1 then

        u := Vectors.householderVector(bh1, cat(
              1,
              fill(0, nx - 1),
              {1}));

        Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);

        bh1 := Vectors.householderReflexion(bh1, u);
        ch1 := Vectors.householderReflexion(ch1, u);

        ll := nx;

        h := true;
        r := 1;

        while r <= nx - 2 and h loop
          if max(Ah1[1:ll - 1, ll]) <= 1e-8 then

            u := cat(
                  1,
                  Vectors.householderVector(Ah1[1:ll - 1, ll], cat(
                    1,
                    fill(0, ll - 2),
                    {1})),
                  fill(0, nx - ll + 1));
            Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);
            ch1 := Vectors.householderReflexion(ch1, u);

            Ah1[1:ll - 2, ll] := fill(0, ll - 2);

            ll := ll - 1;

          else
            h := false;

          end if;
          r := r + 1;
        end while;
      end if;

      Qc2 := cat(
            2,
            Ah1[:, 2:nx],
            matrix(bh1));
      sigma := Modelica.Math.Matrices.singularValues(Qc2);
      eps2 := if eps > 0 then eps else 1000*sigma[1]*Modelica.Constants.eps;
      rankQc2 := 0;

      for i in 1:nx loop
        Modelica.Utilities.Streams.print(" s[" + String(i) + "] = " + String(
          sigma[i]));
        if sigma[i] > eps2 then
          rankQc2 := rankQc2 + 1;
        end if;
      end for;

      Modelica.Utilities.Streams.print("rankQc = " + String(rankQc) +
        "     rankQc2 = " + String(rankQc2) + "     eps2 = " + String(eps2));

      sst := StateSpace(
            A=Ah1,
            B=matrix(bh1),
            C=transpose(matrix(ch1)),
            D=ss.D);

      Ah2 := Ah1[nx - rankQc2 + 1:nx, nx - rankQc2 + 1:nx];
      bh2 := bh1[nx - rankQc2 + 1:nx];
      ch2 := ch1[nx - rankQc2 + 1:nx];
      ssm1 := StateSpace(
            A=Ah2,
            B=matrix(bh2),
            C=transpose(matrix(ch2)),
            D=ss.D);

    end reducedCtrSystem2;

    encapsulated function scaleFactor1
      "Return scale factor for first order block"
      import Modelica;
      input Real n "(s+n)/(s+d)";
      input Real d "(s+n)/(s+d)";
      input Real small=100*Modelica.Constants.eps;
      output Real k "= d/n, if d,n are not zero, otherwise special cases";
    algorithm
      //  k := (if abs(d) > small then abs(d) else 1)/(if abs(n) > small then abs(n) else 1);
      k := if abs(d) > small and abs(n) > small then abs(d)/abs(n) else 1;

      //  k := if abs(n)<=small then 1 else  (if abs(d) > small then abs(d) else 1)/abs(n);

    end scaleFactor1;

    function scaleFactor2 "Return scale factor for second order block"
      import Modelica;
      input Real n1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
      input Real n2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
      input Real d1 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
      input Real d2 "(s^2 + n1*s + n2)/(s^2 + d1*s + d2)";
      input Real small=100*Modelica.Constants.eps;
      output Real k "= d2/n2, if d2,n2 are not zero, otherwise special cases";
    algorithm
      //  k := (if abs(d2) > small then abs(d2) else (if abs(d1) > small then abs(
      //    d1) else 1))/(if abs(n2) > small then abs(n2) else (if abs(n1) > small then
      //          abs(n1) else 1));

      //  if abs(d2) > small and abs(n2) > small then
      //    k := d2/n2;
      //  elseif abs(d2) < small and abs(n2) < small and abs(d1) > small and abs(n1) > small then
      //    k := d1/n1;
      //  else
      //    k := 1;
      //  end if;

      k := if abs(d2) > small and abs(n2) > small then d2/n2 else 1;

    end scaleFactor2;

    encapsulated function staircaseQR
      "Staircase algorithm to put a state space system to controller Hessenberg form"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss "State space system";

      output Boolean isControllable;
      output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "controllable state space system";
      output Real PP[:, :];
    protected
      Real A[:, :];
      Real B[size(ss.B, 1), size(ss.B, 2)];
      Real C[:, :];
      Real Q[:, :];
      Real Q2[:, :];
      Real R[:, :];

      Real P[:, :];
      Real tau[:];

      Integer nx=size(ss.A, 1);
      Integer nu=size(ss.B, 2);
      Integer n1;
      Boolean stop;
      Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);
      Real eps=normA*1e-10;
      Integer stairStepinSys;
      Integer info;
      Integer nn;
      Integer stairStep;
      Integer rankR;

    algorithm
      if nu > 0 then
        if nx > 1 then

          //#####  first step of staircase
          // transform b->Q'b = {*,0,...,0} and c->cQ, A->Q'AQ

          (Q,R,tau,Q2) := Matrices.QR(ss.B);
          B := [R; zeros(nx - nu, nu)];
          // should be the same as transopse(Q2)*ss.B

          A := transpose(Q2)*ss.A;
          A := A*Q2;
          C := ss.C*Q2;
          PP := transpose(Q2);

          stairStep := 0;
          rankR := 0;
          // for i in 1:size(R, 1) loop
          //   if abs(R[i, i]) > eps then
          //     rankR := rankR + 1;
          //   end if;
          // end for;

          //  !!!! rank has to be determined. In the case of ill conditioned systems svd should be used
          for i in 1:min(size(R, 1), size(R, 2)) loop
            if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
              rankR := rankR + 1;
            end if;
          end for;
          stairStep := stairStep + rankR;
          n1 := nx - stairStep;
          stop := false;

          // #######  buildig rest of staircase
          while not stop loop

            (Q,R,tau,Q2) := Matrices.QR(A[stairStep + 1:nx, stairStep - rankR
               + 1:stairStep]);
            P := [identity(nx - n1), zeros(nx - n1, n1); zeros(n1, nx - n1), Q2];
            PP := transpose(P)*PP;
            A := [A[1:stairStep, 1:stairStep], A[1:stairStep, stairStep + 1:nx]
              *Q2; transpose(Q2)*A[stairStep + 1:nx, 1:stairStep], transpose(Q2)
              *A[stairStep + 1:nx, stairStep + 1:nx]*Q2];
            //=transpose(P)*A*P = [A11, A12*Q2; transpose(Q2)*A21, transpose(Q2)*A22*Q2]
            C[:, nx - n1 + 1:nx] := C[:, nx - n1 + 1:nx]*Q2;

            rankR := 0;
            // for i in 1:size(R, 1) loop
            //   if abs(R[i, i]) > eps then
            //     rankR := rankR + 1;
            //   end if;
            // end for;

            //  !!!! rank has to be determined. In the case of ill conditioned systems svd should be used
            for i in 1:min(size(R, 1), size(R, 2)) loop
              if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
                rankR := rankR + 1;
              end if;
            end for;
            stairStep := stairStep + rankR;
            n1 := if rankR < 1 then -1 else n1 - rankR;
            stop := n1 <= 0;

          end while;
        else
          stairStep := if Modelica.Math.Matrices.isEqual(ss.B, zeros(size(ss.B,
            1), size(ss.B, 2))) then 0 else 1;
          A := ss.A;
          B := ss.B;
          C := ss.C;
        end if;

        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=A,
              B=B,
              C=C,
              D=ss.D,
              r=stairStep);

        isControllable := stairStep == nx;
      else
        // no inputs, nu==0
        isControllable := false;
        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=ss.A,
              B=ss.B,
              C=ss.C,
              D=ss.D,
              r=0);
        P := identity(nu);
      end if;

      annotation (Documentation(info="<html>
This algorithm usues QR factorization to generate staircase form i.e. block upper Hessenberg form of the pair (A,B). Due to the well known problem to determine
numerically reliable the rank of a matrix, this algorithm should only be used to well conditioned systems. The best way for rank decision would be singular value decomposition, that is used in staicasSVD.
</html>"));
    end staircaseQR;

    encapsulated function staircaseSVD
      "Staircase algorithm based on singular value decomposition"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output Boolean isControllable;
      output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1), size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "Upper block Hessenberg form state space system";
      output Real P[:, :];

    protected
      Real A[size(ss.A, 1), size(ss.A, 2)];
      Real B[size(ss.B, 1), size(ss.B, 2)];
      Real C[size(ss.C, 1), size(ss.C, 2)];
      Real U[:, :];
      Real VT[:, :];

      Real Q[:, :];
      Real Pi[:, :];
      Real tau[:];
      Real sigma[:];

      Integer nx=size(ss.A, 1);
      Integer nu=size(ss.B, 2);
      Integer ni;
      Boolean stop;
      Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);
      Real eps=normA*1e-10;
      Integer stairStepinSys;
      Integer info;
      Integer nn;
      Integer stairStep;
      Integer rankS;
      Integer maxA=max(size(ss.A));
      //   Real evc[:,2];
      //   Real evnc[:,2];

    algorithm
      if nu > 0 then
        if nx > 1 then
          //#####  first step of staircase
          // transform b->Q'b = {*,0,...,0} and c->cQ, A->Q'AQ
          (sigma,U,VT) := Modelica.Math.Matrices.singularValues(ss.B);

          rankS := 0;
          for i in 1:size(sigma, 1) loop
            if sigma[i] > maxA*sigma[1]*Modelica.Constants.eps then
              rankS := rankS + 1;
            end if;
          end for;

          B := [diagonal(sigma[1:rankS]), zeros(rankS, nu - rankS); zeros(nx -
            rankS, nu)];
          P := transpose(U);
          Q := transpose(VT);
          A := transpose(U)*ss.A*U;
          //    C := ss.C*U;
          (sigma,U,VT) := Modelica.Math.Matrices.singularValues(A[rankS + 1:nx,
            1:rankS]);

          stairStep := rankS;
          rankS := 0;
          if size(sigma, 1) > 1 then
            for i in 1:size(sigma, 1) loop
              if sigma[i] > maxA*sigma[1]*Modelica.Constants.eps then
                rankS := rankS + 1;
              end if;
            end for;
          else
            rankS := if size(sigma, 1) > 0 then if sigma[1] > maxA*Modelica.Constants.eps
               then 1 else 0 else 0;
          end if;

          stairStep := stairStep + rankS;
          Pi := [VT, zeros(size(VT, 1), size(U, 1)); zeros(size(U, 2), size(VT,
            2)), transpose(U)];
          B := Pi*B;
          P := Pi*P;

          // P could be ambigious according to the sign
          if transpose(P)*B[:, 1]*ss.B[:, 1] < 0 then
            Pi := -Pi;
            P := -P;
          end if;

          A := Pi*A*transpose(Pi);

          // should be made better because of many zeros in B

          while stairStep < nx and rankS > 0 and not
              Modelica.Math.Matrices.isEqual(
                  A[stairStep + 1:nx, stairStep - rankS + 1:stairStep],
                  zeros(stairStep - rankS, rankS),
                  eps) loop

            (sigma,U,VT) := Modelica.Math.Matrices.singularValues(A[stairStep
               + 1:nx, stairStep - rankS + 1:stairStep]);

            Pi := [identity(stairStep - rankS), zeros(stairStep - rankS, nx -
              stairStep + rankS); zeros(nx - stairStep + rankS, stairStep -
              rankS), [VT, zeros(rankS, nx - stairStep); zeros(nx - stairStep,
              rankS), transpose(U)]];
            P := Pi*P;
            A := Pi*A*transpose(Pi);

            //new implenmentation advisable because of many zeros in Pi
            //    C := C*transpose(Pi);
            rankS := 0;
            if size(sigma, 1) > 1 then
              for i in 1:size(sigma, 1) loop
                if sigma[i] > maxA*sigma[1]*Modelica.Constants.eps then
                  rankS := rankS + 1;
                end if;
              end for;
            else
              rankS := if size(sigma, 1) > 0 then if sigma[1] > maxA*Modelica.Constants.eps
                 then 1 else 0 else 0;
            end if;
            stairStep := stairStep + rankS;
          end while;

          B := P*ss.B;
          C := ss.C*transpose(P);

        else
          stairStep := if Modelica.Math.Matrices.isEqual(ss.B, zeros(size(ss.B,
            1), size(ss.B, 2))) then 0 else 1;
          A := ss.A;
          B := ss.B;
          C := ss.C;
        end if;

        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=A,
              B=B,
              C=C,
              D=ss.D,
              r=stairStep);
        isControllable := stairStep == nx;

      else
        // no inputs, nu==0
        isControllable := false;
        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=ss.A,
              B=ss.B,
              C=ss.C,
              D=ss.D,
              r=0);
        P := identity(nu);
      end if;

    end staircaseSVD;

    encapsulated function trianUpperHess
      "Triangulize an upper Hessenberg matrix by repeatedly applicated householder reflexion"
      import Modelica;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input Real H[:, :] "Upper Hessenberg matrix";
      input Real b[size(H, 1)];

      output Real Ht[size(H, 1), size(H, 2)];
      output Real bt[size(b, 1)];

    protected
      Integer q=size(H, 1);
      Real u[:] "Householder vector";
      Integer ll;

    algorithm
      Ht := H;
      bt := b;

      for ll in 1:q - 1 loop
        u := cat(
              1,
              zeros(ll - 1),
              cat(
                1,
                Vectors.householderVector(vector(Ht[ll:ll + 1, ll]), {1,0}),
                zeros(q - ll - 1)));
        Ht := Matrices.householderReflexion(Ht, u);
        bt := Vectors.householderReflexion(bt, u);
      end for;

      annotation (Documentation(info="<html>
<p>
This function computes an triangular matrix from an upper Hessenberg matrix by stepwise annihilation of the subdiagonal elements.
</p>
<blockquote><pre>
<b>A</b> -> <b>QA</b> = <b>T</b>
</pre></blockquote>
<p>
It is assumend that the original matrix has upper hessenberg form.
Additionally the vector b is transformed in the same way
</p>
<blockquote><pre>
<b>b</b> -> <b>Qb</b> = <b>q</b>
</pre></blockquote>
<p>
The function is primarily used to calculate the transfer function gain from a SISO state space system in observer Hessenberg form
</p>
<blockquote><pre>
    ( *   *   ...   ...    * )          ( * )
    ( *   *   ...   ...    * )          ( . )
<b>A</b> = ( 0   *   ...   ...    * ),    <b>b</b> =  ( . ),   <b>c</b> = ( 0, ..., 0, * )
    ( .   .    .     .     . )          ( * )
    ( 0  ...   0     *     * )          ( * )
</pre></blockquote>
<p>
If <b>A</b> is upper Hessenberg and <b>T</b> = <b>Q</b>*<b>A</b> is triangular then obviously <b>H</b>(s) = <b>Q</b>*(s*<b>I</b> -<b>A</b>) = s*<b>I</b> - <b>T</b>.
</p>
<p>
Further on, if <b>T</b> is triangular then also <b>H</b> = s<b>I</b> - <b>T</b> is and the element l_nn of <b>L</b> = inv(<b>H</b>) is given by 1/h_nn.
The frequency response G(s0)for a given s0 that is neither zero nor pole of the system can be calculated by
</p>
<blockquote><pre>
G(s0)  = <b>c</b>*(s0*<b>I</b> -<b>A</b>)<sup>-1</sup>*<b>b</b> = <b>c</b>*(s0*<b>I</b> -<b>A</b>)<sup>-1</sup> *<b>Q</b>*<b>Q</b><sup>-1</sup>*<b>b</b> = <b>c</b>*(<b>Q</b><sup>-1</sup>*(s0*<b>I</b> -<b>A</b>))<sup>-1</sup>*<b>Q</b><sup>-1</sup>*<b>b</b> = <b>c</b>*<b>H</b><sup>-1</sup>(s0)*<b>q</b>
</pre></blockquote>
<p>
and because only the n'th element of <b>c</b> is different to zero the gain k is given by
</p>
<blockquote><pre>
    q_nn*c_nn     product(s0 - poles_i)
k = ---------- * ----------------------
       h_nn       product(s0 - zeros_i)
</pre></blockquote>
</html>"));
    end trianUpperHess;

    encapsulated function transposeStateSpace
      "Return the transposed state space system"

      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss "State space system";

      output StateSpace sst=StateSpace(
              A=transpose(ss.A),
              B=transpose(ss.C),
              C=transpose(ss.B),
              D=transpose(ss.D),
              uNames=ss.yNames,
              yNames=ss.uNames);
    algorithm

      annotation (Documentation(info="<html>
</html>"));
    end transposeStateSpace;

    encapsulated function reduceRosenbrock
      "Algorithm to compress the generalized system matrix [A, B; C, D] to calculate the invariant zeros of a system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      import Modelica.Utilities.Streams.print;

      input Real A[:, :];
      input Real B[:, :];
      input Real C[:, :];
      input Real D[:, :];

      output Real Ar[:, :];
      output Real Br[:, :];
      output Real Cr[:, :];
      output Real Dr[:, :];
      output Integer n "= dimension of Ar: Ar[n,n]";
      output Integer m "= second dimension of Br: Br{n,m]";
      output Integer p "= first dimension of Cr: Cr[p,n]";

    protected
      Real A2[:, :];
      Real B2[:, :];
      Real C2[:, :];
      Real CC[:, :];
      Real Co[:, :];
      Real Cu[:, :];
      Real D2[:, :];
      Real DD[:, :];
      Real Mr[:, :];
      Real Vf[:, :];
      Real V[:, :];
      Real V2[:, :];
      Real R[:, :];
      Real tau[:];

      Integer nx=size(A, 1);
      Integer nu=size(B, 2);
      Integer ny=size(C, 1);

      Integer nue;
      Integer delta;
      Integer rho;
      Integer mue;
      Integer sigma;
      Integer j;
      Boolean stop;
      Boolean stop1 "reduction finished";
      Boolean stop2 "system has no zeros";
      Integer rankR;
      // Real normA=Modelica.Math.Matrices.norm(A=A, p=1);
      // Real eps=normA*1e-12;
      Real normABCD=Modelica.Math.Matrices.norm([A, B; C, D], p=1);
      Real eps=normABCD*Modelica.Constants.eps*1000;
    algorithm
      if nx > 0 then

        A2 := A;
        B2 := B;
        C2 := C;
        D2 := D;
        stop := false;
        stop1 := false;
        stop2 := false;
        nue := nx;
        delta := 0;
        mue := ny;
        sigma := ny;
        j := 1;

        while not stop loop
          (V,R,tau,V2) := Matrices.QR(D2);

          rankR := 0;
          //  !!!! rank has to be determined. In the case of ill conditioned systems svd should be used
          for i in 1:min(size(R, 1), size(R, 2)) loop
            if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
              rankR := rankR + 1;
            end if;
          end for;

          //rankR:=Modelica.Math.Matrices.rank(R);

          DD := R[1:rankR, :];

          CC := transpose(V2)*C2;

          sigma := rankR;
          stop1 := size(CC, 1) == rankR;

          if not stop1 then
            Cu := CC[sigma + 1:end, :];
            Co := CC[1:sigma, :];

            (V,R,tau,V2) := Matrices.QR(Matrices.fliplr(transpose(Cu)));
            Vf := Matrices.fliplr(V2);

            rankR := 0;
            //  !!!! rank determination
            for i in 1:min(size(R, 1), size(R, 2)) loop
              if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
                rankR := rankR + 1;
              end if;
            end for;
            //rankR:=Modelica.Math.Matrices.rank(R);

            rho := rankR;
            stop1 := rho == 0;
            stop2 := size(Cu, 2) == rankR;

            if not stop1 and not stop2 then
              nue := size(Cu, 2) - rankR;
              mue := rho + sigma;
              delta := delta + rho;

              if sigma == 0 then
                Mr := [transpose(Vf)*A2*Vf, transpose(Vf)*B2];
              else
                Mr := [transpose(Vf)*A2*Vf, transpose(Vf)*B2; Co*Vf, DD];
              end if;

              A2 := Mr[1:nue, 1:nue];
              B2 := Mr[1:nue, nue + rho + 1:nue + rho + nu];
              C2 := Mr[nue + 1:nue + mue, 1:nue];
              D2 := Mr[nue + 1:nue + mue, nue + rho + 1:nue + rho + nu];

              j := j + 1;
            end if;
            //not stop1 or not stop2

          end if;
          //if not stop1

          stop := stop1 or stop2 or j > 3*nx;

        end while;

        if stop1 then
          Ar := A2;
          Br := B2;
          Cr := C2;
          Dr := D2;
          n := nue;
          p := sigma;
          m := nu;
        else
          n := 0;
          p := 0;
          m := 0;
          Ar := fill(
                0,
                0,
                0);
          Br := fill(
                0,
                0,
                0);
          Cr := fill(
                0,
                0,
                0);
          Dr := fill(
                0,
                0,
                0);
        end if;

      else
        n := 0;
        p := 0;
        m := 0;
        A2 := A;
        B2 := B;
        C2 := C;
        D2 := D;
      end if;

      annotation (Documentation(info="<html>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Emami-Naeini, A. and Van Dooren, P. (1982):</dt>
<dd> <b>Computation of Zeros of Linear Multivariable Systems</b>.
     Automatica, 18, pp. 415-430.<br>&nbsp;</dd>
</dl>
</html>


"));
    end reduceRosenbrock;

    encapsulated function reducedCtrSystem
      "Calculate the controllable part of a SISO system"
      import Modelica.Utilities.Streams.print;
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss "State space system";
      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), 1],
        redeclare Real C[1, size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "controllable state space system";

    protected
      Integer nx=size(ss.A, 1);
      Real Ah1[nx, size(ss.A, 2)];
      Real bh1[nx];
      Real ch1[nx];
      Real u[:] "householder vector";
      Integer ll;
      Integer r=1;
      Real maxa;
      Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      Ah1 := ss.A;
      bh1 := ss.B[:, 1];
      ch1 := ss.C[1, :];

      if Modelica.Math.Vectors.length(bh1) > 0 then
        r := 1;
        if nx > 1 then

          u := Vectors.householderVector(bh1, cat(
                1,
                fill(0, nx - 1),
                {1}));

          Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);
          bh1 := Vectors.householderReflexion_en(bh1, u);
          ch1 := Vectors.householderReflexion(ch1, u);
          bh1[1:nx - 1] := fill(0, nx - 1);
          ll := nx;
          maxa := max(abs(Ah1[1:ll - 1, ll]));

          while r <= nx - 1 and maxa > normA*tol loop
            u := cat(
                  1,
                  Vectors.householderVector(Ah1[1:ll - 1, ll], cat(
                    1,
                    fill(0, ll - 2),
                    {1})),
                  fill(0, nx - ll + 1));

            Ah1 := Matrices.Internal.hohoTrafoLowerHess(
                  Ah1,
                  u,
                  r);

            ch1 := Vectors.householderReflexion(ch1, u);
            ll := ll - 1;
            maxa := max(abs(Ah1[1:ll - 1, ll]));

            r := r + 1;
          end while;
        end if;
        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=Ah1,
              B=matrix(bh1),
              C=transpose(matrix(ch1)),
              D=ss.D,
              r=r);
        ssm1.r := r;
      end if;
    end reducedCtrSystem;

    encapsulated function reducedCtrSystemX
      "Calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

      input StateSpace ss "State space system";

      output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), 1],
        redeclare Real C[1, size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1), size(ss.D, 2)])
        "Controllable state space system";

    protected
      Integer nx=size(ss.A, 1);
      Real Ah1[nx, size(ss.A, 2)];
      Real bh1[nx];
      Real ch1[nx];
      Real u[:] "householder vector";
      Integer ll;
      Integer r=1;
      Real maxa;
      Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);

    algorithm
      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1,
          "A SISO-system is expected as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1,
          " A SISO-system is expected as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;

      Ah1 := ss.A;
      bh1 := ss.B[:, 1];
      ch1 := ss.C[1, :];

      if Modelica.Math.Vectors.length(bh1) > 0 then
        r := 1;
        if nx > 1 then

          u := Vectors.householderVector(bh1, cat(
                1,
                fill(0, nx - 1),
                {1}));

          Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);
          bh1 := Vectors.householderReflexion_en(bh1, u);
          ch1 := Vectors.householderReflexion(ch1, u);
          bh1[1:nx - 1] := fill(0, nx - 1);

          ll := nx;
          maxa := max(abs(Ah1[1:ll - 1, ll]));

          while r <= nx - 1 and maxa > nx*normA*1e-5 loop
            u := cat(
                  1,
                  Vectors.householderVector(Ah1[1:ll - 1, ll], cat(
                    1,
                    fill(0, ll - 2),
                    {1})),
                  fill(0, nx - ll + 1));

            Ah1 := Matrices.Internal.hohoTrafoLowerHess(
                  Ah1,
                  u,
                  r);

            ch1 := Vectors.householderReflexion(ch1, u);
            ll := ll - 1;
            maxa := max(abs(Ah1[1:ll - 1, ll]));

            r := r + 1;
          end while;

        end if;

        ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
              A=Ah1,
              B=matrix(bh1),
              C=transpose(matrix(ch1)),
              D=ss.D,
              r=r);
        ssm1.r := r;
      end if;
    end reducedCtrSystemX;

    function read_dslin "Read a StateSpace data record from mat-file"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

      input String fileName = "dslin" "Name of the result file";
        //annotation (Dialog(
        //    loadSelector(filter="MAT files (*.mat);; All files (*.*)",
        //    caption="State space system data file")));

    protected
      String fileName2 = fileName + ".mat" "Name of the result file with extension";
      Integer xuy[3] = ReadSystemDimension(fileName2, "ABCD");
      Integer nx = xuy[1];
      Integer nu = xuy[2];
      Integer ny = xuy[3];
      Real ABCD[nx + ny, nx + nu] = Streams.readRealMatrix(fileName2, "ABCD", nx + ny, nx + nu);
      String xuyName[nx + nu + ny] = readStringMatrix(fileName2, "xuyName",nx + nu + ny);
    public
      output StateSpace result(
        redeclare Real A[nx, nx],
        redeclare Real B[nx, nu],
        redeclare Real C[ny, nx],
        redeclare Real D[ny, nu]) "Outputs model linearized at initial point";

    algorithm
      result.A := ABCD[1:nx, 1:nx];
      result.B := ABCD[1:nx, nx + 1:nx + nu];
      result.C := ABCD[nx + 1:nx + ny, 1:nx];
      result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
      result.uNames := xuyName[nx + 1:nx + nu];
      result.yNames := xuyName[nx + nu + 1:nx + nu + ny];
      result.xNames := xuyName[1:nx];

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = StateSpace.Import.<b>fromModel</b>(modelName, T_linearize, fileName)
</pre></blockquote>

<h4>Description</h4>
<p>Generate a StateSpace data record by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation. The result of linearization is transformed into a StateSpace record.
</p>

<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
  Real T_linearize = 5;

<b>algorithm</b>
  ss = Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName, T_linearize);

// ss.A = [ 0.0,   1.0,    0.0,            0.0,      0.0,     0.0;
            0.0,   0.0,          -2.26,    0.08,     1.95,   -0.45;
            0.0,   0.0,           0.0,            1.0,      0.0,     0.0;
            0.0,   0.0,          -3.09,   -1.38,     7.70,   -3.01;
            0.0,   0.0,           0.0,            0.0,      0.0,     1.0;
            0.0,   0.0,          -6.47,    1.637,   -2.90,    1.29],

// ss.B=[0.0; 0.13; 0.0; -0.014; 0.0; -0.1],
// ss.C=identity(6),
// ss.D=[0; 0; 0; 0; 0; 0]
</pre></blockquote>
</html>"));
    end read_dslin;

    encapsulated function absComplexVector
      "Return the absolute values of all elements of a complex vector that is defined by a vr[:] vector (real part) and a vi[:] vector (imaginary part)"
      import Modelica;
      input Real vr[:] "Real part of complex vector";
      input Real vi[size(vr,1)] "Imaginary part of complex vector";
      output Real v_abs[size(vr,1)]
        "Absolute values of the elements of the complex vector";
    protected
      Real r_abs;
      Real i_abs;
    algorithm
       for i in 1:size(vr,1) loop
          r_abs :=abs(vr[i]);
          i_abs :=abs(vi[i]);
          v_abs[i] :=if r_abs == 0 and i_abs == 0 then 0 else if r_abs > i_abs
           then r_abs*sqrt(1 + (i_abs/r_abs)^2) else i_abs*sqrt(1 + (r_abs/i_abs)^2);
       end for;
    end absComplexVector;

    encapsulated function polesAndZeros_Old
      "Plot poles (i.e. eigenvalues) and/or invariant zeros of a state space system (previous version of polesAndZeros that is kept, just in case)"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "Linear system in state space form"
        annotation (Dialog);
      input Boolean poles=true
        "= true, to plot the poles (i.e. the eigenvalues) of ss"
        annotation (choices(checkBox=true));
      input Boolean zeros=true "= true, to plot the (invariant) zeros of ss "
        annotation (choices(checkBox=true));

      input Boolean print=false
        "= true, to print the selection to the output window"
        annotation (choices(checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=if poles and zeros then
            Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros()
             else if poles then
            Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
            heading="Poles (x)") else
            Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros(
            heading="Invariant zeros (o)"));
    protected
      Integer nx=size(ss.A, 1);
      Real eval[nx, 2];
      Real invZerosRe[:];
      Real invZerosIm[:];
      Complex invZeros[:];
      Complex eig[size(ss.A,1)];
      Plot.Records.Curve curves[2];
      Integer i;
      Plot.Records.Diagram diagram2;
    algorithm
      // Determine eigen values
      if poles then
        eval := Modelica.Math.Matrices.eigenValues(ss.A);
        //eval :=Modelica_LinearSystems2.Math.Matrices.eigenValuesAsRealMatrix(ss.A);
      end if;

      if zeros then
        invZeros := StateSpace.Analysis.invariantZeros(ss);
        invZerosRe := fill(0, size(invZeros, 1));
        invZerosIm := fill(0, size(invZeros, 1));
        for i in 1:size(invZeros, 1) loop
          invZerosRe[i] := invZeros[i].re;
          invZerosIm[i] := invZeros[i].im;
        end for;
      end if;

      i := 0;
      if poles then
        i := i + 1;
        curves[i] := Plot.Records.Curve(
              x=eval[:, 1],
              y=eval[:, 2],
              legend="poles",
              autoLine=false,
              linePattern=Plot.Types.LinePattern.None,
              lineSymbol=Plot.Types.PointSymbol.Cross);
      end if;

      if zeros then
        i := i + 1;
        curves[i] := Plot.Records.Curve(
              x=invZerosRe,
              y=invZerosIm,
              legend="zeros",
              autoLine=false,
              linePattern=Plot.Types.LinePattern.None,
              lineSymbol=Plot.Types.PointSymbol.Circle);
      end if;

      diagram2 := defaultDiagram;
      diagram2.curve := curves[1:i];
      Plot.diagram(diagram2, device);

      if print then
         if poles then
            for i in 1:nx loop
              eig[i].re := eval[i,1];
              eig[i].im := eval[i,2];
            end for;
          Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.printHTML(
            eig, heading="Eigenvalues", name="eigenvalue");
         end if;

         if zeros then
          Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.printHTML(
            invZeros, heading="Invariant zeros", name="invariant zero");
         end if;
      end if;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>polesAndZeros</b>(ss);
   or
StateSpace.Plot.<b>polesAndZeros</b>(
  ss,
  poles=true,
  zeros=true,
  plot=true,
  defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros</a>(),
  device=<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>());
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a pole-zero-map of the poles and transmission zeros of a state space system.
The poles are the eigenvalues of the system matrix (eigenvalues(ss.A)). The Boolean inputs
\"poles\" and \"zeros\" define what to plot. If Boolean input \"plot = true\", the pole-zero-map
is plotted. If false, only the diagram is generated and returned as output argument.
The records \"defaultDiagram\" and \"device\" allow to set various layout options and the
size and location of the diagram on the screen.
</p>

<h4>Example</h4>
<p>
The example <a href=\"modelica://Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros\">
Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros</a>
is defined as
</p>
<pre>
  Plot.polesAndZeros(ss = Modelica_LinearSystems2.StateSpace(
    A=[-3, 2,-3,  4, 5,6;
        0, 6, 7,  8, 9,4;
        0, 2, 3,  0,78,6;
        0, 1, 2,  2, 3,3;
        0,13,34,  0, 0,1;
        0, 0, 0,-17, 0,0],
    B=[1,0;
       0,1;
       1,0;
       0,1;
       1,0;
       0,1],
    C=[0,0,1,0,1,0;
       0,1,0,0,1,1],
    D=[0,0;
       0,0]));
</pre>

<p>
and results in
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/StateSpace/polesAndZerosSS.png\"/>
</blockquote>
</html>"));
    end polesAndZeros_Old;

    encapsulated function bodeSISO_Old
      "Plot bode plot of the corresponding transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica.Utilities.Streams.print;

      input StateSpace ss "State space system";
      input Integer iu=1 "Index of input";
      input Integer iy=1 "Index of output";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange=true
        "True, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min=0.1
        "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max=10
        "Maximum frequency value, if autoRange = false";

      input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
      input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
          defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

      input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));
      input Boolean onFile=false
        "= true, if frequency response is stored on file as matrix [f,A,phi]" annotation(choices(checkBox=true));
      input String fileName="frequencyResponse.mat"
        "If onFile=true, file on which the frequency response will be stored"  annotation(Dialog(enable=onFile));
      input String matrixName=if Hz and not dB then "fHz_A_phiDeg" elseif
                                 Hz and dB then "fHz_AdB_phiDeg" elseif
                                 not Hz and dB then "f_AdB_phiDeg" else "f_A_phiDeg"
        "If onFile=true, Name of matrix on file" annotation(Dialog(enable=onFile));

    protected
      ZerosAndPoles zp "ZP-Transfer functions to be plotted";
      StateSpace ss_siso(
        redeclare Real A[size(ss.A, 1), size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1), 1],
        redeclare Real C[1, size(ss.C, 2)],
        redeclare Real D[1, 1]);

    algorithm
      // Check that system has inputs and outputs
      if size(ss.B, 2) == 0 then
        Modelica.Utilities.Streams.print(
          "\n... Not possible to plot transfer function because system has no inputs."
           + "\n... Call of Plot.bodeSISO is ignored.\n");
        return;
      elseif size(ss.C, 1) == 0 then
        Modelica.Utilities.Streams.print(
          "\n... Not possible to plot transfer function because system has no outputs."
           + "\n... Call of Plot.bodeSISO is ignored.\n");
        return;
      end if;

      assert(iu <= size(ss.B, 2) and iu > 0, "index for input is " + String(iu)
         + " which is not in [1, " + String(size(ss.B, 2)) + "].");
      assert(iy <= size(ss.C, 1) and iy > 0, "index for output is " + String(iy)
         + " which is not in [1, " + String(size(ss.C, 1)) + "].");
      ss_siso := StateSpace(
            A=ss.A,
            B=matrix(ss.B[:, iu]),
            C=transpose(matrix(ss.C[iy, :])),
            D=matrix(ss.D[iy, iu]));

      zp := StateSpace.Conversion.toZerosAndPoles(
               StateSpace.Transformation.toBalancedForm(ss_siso), tol);

      ZerosAndPoles.Plot.bode(
            zp,
            nPoints,
            autoRange,
            f_min,
            f_max,
            Hz=Hz,
            magnitude=magnitude,
            dB=dB,
            phase=phase,
            onFile=onFile,
            fileName=fileName,
            matrixName=matrixName,
            defaultDiagram=defaultDiagram,
            device=device);

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>bodeSISO</b>(ss)
   or
StateSpace.Plot.<b>bodeSISO</b>(
  ss,
  iu,
  iy,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the bode-diagram of a transfer function corresponding
to the behavior of the state space system from iu'th element of the input
vector <b>u</b> to the iy'th element of the output vector <b>y</b>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0],
    B=[0.0,1.0; 1.0,1.0; -1.0,0.0],
    C=[0.0,1.0,1.0; 1.0,1.0,1.0],
    D=[1.0,0.0; 0.0,1.0])

  Integer iu=1;
  Integer iy=1;

<b>algorithm</b>
   Modelica_LinearSystems2.StateSpace.Plot.plotBodeSISO(ss, iu, iy)
//  gives:
</pre></blockquote>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</p>
</html>"));
    end bodeSISO_Old;

    encapsulated function bodeMIMO_old
      "Plot bode plot of all transfer functions, corresponding to the state space system"
      import Modelica.Utilities.Streams.print;
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "State space system";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange[:, :]=fill(
          true,
          size(ss.C, 1),
          size(ss.B, 2)) "True, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min[:, :]=fill(
          0.1,
          size(ss.C, 1),
          size(ss.B, 2)) "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max[:, :]=fill(
          10,
          size(ss.C, 1),
          size(ss.B, 2)) "Maximum frequency value, if autoRange = false";

      input Boolean magnitude=true "= true, to plot magnitude" annotation(choices(checkBox=true));
      input Boolean phase=true "= true, to plot phase" annotation(choices(checkBox=true));

      input Real tol=1e-10
        "Tolerance of reduction procedure, default tol = 1e-10";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

      input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)" annotation(choices(checkBox=true));
      input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))" annotation(choices(checkBox=true),Dialog(enable=magnitude));
    protected
      ZerosAndPoles zp[size(ss.C, 1), size(ss.B, 2)]
        "ZerosAndPoles object to be plotted";
      Plot.Records.Diagram diagram2=defaultDiagram;
      String yNames[size(ss.C, 1)];
      String uNames[size(ss.B, 2)];

    algorithm
      // Check that system has inputs and outputs
      if size(ss.B, 2) == 0 then
        Modelica.Utilities.Streams.print("\n... Not possible to plot transfer function because system has no inputs."
           + "\n... Call of Plot.bodeMIMO is ignored.\n");
        return;
      elseif size(ss.C, 1) == 0 then
        Modelica.Utilities.Streams.print("\n... Not possible to plot transfer function because system has no outputs."
           + "\n... Call of Plot.bodeMIMO is ignored.\n");
        return;
      end if;

      // generate headings
      for i1 in 1:size(ss.B, 2) loop
        uNames[i1] := if ss.uNames[i1] == "" then "u" + String(i1) else ss.uNames[
          i1];
      end for;
      for i1 in 1:size(ss.C, 1) loop
        yNames[i1] := if ss.yNames[i1] == "" then "y" + String(i1) else ss.yNames[
          i1];
      end for;

      zp := StateSpace.Conversion.toZerosAndPolesMIMO(
            StateSpace.Transformation.toBalancedForm(ss), tol);

      for i1 in 1:size(ss.C, 1) loop
        for i2 in 1:size(ss.B, 2) loop
          diagram2.heading := defaultDiagram.heading + "  " + uNames[i2] + " -> " +
            yNames[i1];
          ZerosAndPoles.Plot.bode(
            zp[i1, i2],
            nPoints,
            autoRange[i1, i2],
            f_min[i1, i2],
            f_max[i1, i2],
            magnitude=magnitude,
            phase=phase,
            Hz=Hz,
            dB=dB,
            defaultDiagram=diagram2,
            device=device);
        end for;
      end for;

      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
StateSpace.Plot.<b>bodeMIMO</b>(ss)
   or
StateSpace.Plot.<b>bodeMIMO</b>(
  ss,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0],
    B=[0.0,1.0; 1.0,1.0; -1.0,0.0],
    C=[0.0,1.0,1.0],
    D=[1.0,0.0])

<b>algorithm</b>
   Modelica_LinearSystems2.StateSpace.Plot.plotBodeMIMO(ss)
//  gives:
</pre></blockquote>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase.png\">
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagnitude2.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhase2.png\">
</p>
</html>"));
    end bodeMIMO_old;
  end Internal;

  annotation (defaultComponentName="stateSpace", Documentation(info="<html>
<p>
This record defines a linear time invariant differential
equation system in state space form:
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>
</pre></blockquote>
<p>
with
</p>
<ul>
<li> <b>u</b> ... the input vector,</li>
<li> <b>y</b> ... the output vector,</li>
<li> <b>x</b> ... the state vector,</li>
<li> <b>A</b>, <b>B</b>, <b>C</b>, <b>D</b> - matrices of appropriate dimensions.</li>
</ul>
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
          lineColor={255,255,255},
          extent={{-90,-50},{90,50}},
          textString="SS")}));
end StateSpace;
