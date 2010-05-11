within Modelica_LinearSystems2;
record StateSpace
  "Continuous state space description of a linear, time invariant differential equation system (data + operations)"

  extends Modelica.Icons.Record;

  Real A[:,size(A, 1)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real B[size(A, 1),:]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real C[:,size(A, 1)]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));

   String yNames[size(C, 1)]=fill("", size(C, 1)) "Names of the output signals"
                                                                                annotation(Dialog(group="Signal names"));
   String xNames[size(A, 1)]=fill("", size(A, 1)) "Names of the states"  annotation(Dialog(group="Signal names"));
   String uNames[size(B, 2)]=fill("", size(B, 2)) "Names of the input signals"  annotation(Dialog(group="Signal names"));

encapsulated operator 'constructor'
    "Default constructors for a StateSpace record"
    import Modelica_LinearSystems2;

  function fromABCDMatrices "Default constructor for a StateSpace record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input Real A[:,size(A, 1)];
      input Real B[size(A, 1),:];
      input Real C[:,size(A, 1)];
      input Real D[size(C, 1),size(B, 2)];

      input String uNames[size(B, 2)]=fill("", size(B, 2));
      input String yNames[size(C, 1)]=fill("", size(C, 1));
      input String xNames[size(A, 2)]=fill("", size(A, 2));

      output StateSpace result(
        redeclare Real A[size(A, 1),size(A, 2)],
        redeclare Real B[size(B, 1),size(B, 2)],
        redeclare Real C[size(C, 1),size(C, 2)],
        redeclare Real D[size(D, 1),size(D, 2)],
        redeclare String uNames[size(B, 2)],
        redeclare String yNames[size(C, 1)],
        redeclare String xNames[size(A, 2)]);

  algorithm
      result.A := A;
      result.B := B;
      result.C := C;
      result.D := D;
      result.uNames := uNames;
      result.yNames := yNames;
      result.xNames := xNames;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center>=</td>  <td> 'constructor'.<b>fromABCDMatrices</b>(A, B, C, D)  </td> </tr>

</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function constructs a StateSpace record ss with<br>
<blockquote><pre>
  ss.A = A;
  ss.B = B;
  ss.C = C;
  ss.D = D;
</pre></blockquote>

</p>

<h4><font color=\"#008000\">Example</font></h4>
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

  function fromReal "Generate a StateSpace data record from a Real value"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

    input Real r "Value of Real variable";
    output StateSpace ss(
      redeclare Real A[0,0],
      redeclare Real B[0,1],
      redeclare Real C[1,0],
      redeclare Real D[1,1]) "= r";

  algorithm
    ss.D[1, 1] := r;
    annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center>=</td>  <td> 'constructor'.<b>fromReal</b>(r)  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function constructs a StateSpace record ss from a Real value, i.e. a state space system without a state and an output without dynamics:
<blockquote><pre>
y = r*u
</pre></blockquote>
Therefore, the matrices are defined by
<blockquote><pre>
  ss.A = fill(0,0,0);
  ss.B = fill(0,0,1);
  ss.C = fill(0,1,0);
  ss.D = [r];
</pre></blockquote>
 
</p>
 
 
</html>"));
  end fromReal;

  function fromTransferFunction = 
      Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace annotation (Documentation(info="<html> </html>"));
  function fromZerosAndPoles = 
      Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace annotation (Documentation(info="<html> </html>"));

    annotation (Documentation(info="<html>
This package contains the default constructors for StateSpace record.
</html>"));
end 'constructor';

encapsulated operator '-'
    "Contains operators for subtraction of state space systems"

  function subtract
      "Subtraction of two state space systems connected in parallel (= inputs are the same, outputs of the two systems are subtracted)"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "State space system 1";
    input StateSpace ss2 "State Space system 2 is subtracted from system 1";
    output StateSpace result(
      redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1),size(ss1.A, 2) + size(
        ss2.A, 2)],
      redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1),size(ss1.B, 2)],
      redeclare Real C[size(ss1.C, 1),size(ss1.C, 2) + size(ss2.C, 2)],
      redeclare Real D[size(ss1.D, 1),size(ss1.D, 2)]) "= ss1 - ss2";
    protected
    Integer nx1=size(ss1.A, 1);
    Integer nx2=size(ss2.A, 1);
  algorithm
    result.A := [ss1.A,zeros(nx1, nx2); zeros(nx2, nx1),ss2.A];
    result.B := [ss1.B; ss2.B];
    result.C := [ss1.C,-ss2.C];
    result.D := ss1.D - ss2.D;
      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> Modelica_LinearSystems2.StateSpace.'-'.<b>subtract</b>(ss1, ss2)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This operator function computes the subtraction of two state space systems connected in parallel, i.e. the inputs are the same and the outputs of the two systems are subtracted. Therefore, The systems must have the same number of inputs and outputs but not the same number of states. The resulting system has an order of system_order1 + system_order2.<br>
The operator is used by writing just the following command:
<blockquote><pre>
    ss3 := ss1 - ss2;
</pre></blockquote>

</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   StateSpace ss1 = StateSpace(A=[-1, 0; 0, -2], B=[1;2], C=[0, 1], D=[0]);
   StateSpace ss2 = StateSpace(A=[-3, 0; 0, -4], B=[3;4], C=[0, 2], D=[0]);
   
   StateSpace ss3;

<b>algorithm</b>
  ss3 := ss1 - ss2;
// ss.A = [-1, 0, 0, 0; 0, -2, 0, 0; 0, 0, -3, 0; 0, 0, 0, -4],
// ss.B = [1; 2; 3; 4],
// ss.C = [0, 1, 0, -2],
// ss.D = [0],
</pre></blockquote>

</html> "));
  end subtract;

  function negate
      "Unary minus (state space system where the output is multiplied by a gain of -1)"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;
    output StateSpace result(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
      redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)]) "= -ss";
  algorithm
    result.A := ss.A;
    result.B := ss.B;
    result.C := -ss.C;
    result.D := -ss.D;
  end negate;
    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This package contains the <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.'-'.subtract\">'subtract'</a> and the <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.'-'.subtract\">'-'</a> operator for StateSpace records. 

</html>"));
end '-';

encapsulated operator function '+'
    "Parallel connection of two state space systems (= inputs are the same, outputs of the two systems are added)"
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "System 1";
    input StateSpace ss2 "System 2 is added in parallel to system 1";
    output StateSpace result(
      redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1),size(ss1.A, 2) + size(
        ss2.A, 2)],
      redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1),size(ss1.B, 2)],
      redeclare Real C[size(ss1.C, 1),size(ss1.C, 2) + size(ss2.C, 2)],
      redeclare Real D[size(ss1.D, 1),size(ss1.D, 2)]) "= ss1 + ss2";
  protected
    Integer nx1=size(ss1.A, 1);
    Integer nx2=size(ss2.A, 1);
algorithm
    result.A := [ss1.A,zeros(nx1, nx2); zeros(nx2, nx1),ss2.A];
    result.B := [ss1.B; ss2.B];
    result.C := [ss1.C,ss2.C];
    result.D := ss1.D + ss2.D;

end '+';

encapsulated operator function '*'
    "Series connection of two state space systems"
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "System 1";
    input StateSpace ss2 "System 2";
    output StateSpace result(
      redeclare Real A[size(ss1.A, 1) + size(ss2.A, 1),size(ss1.A, 2) + size(
        ss2.A, 2)],
      redeclare Real B[size(ss1.B, 1) + size(ss2.B, 1),size(ss2.B, 2)],
      redeclare Real C[size(ss1.C, 1),size(ss1.C, 2) + size(ss2.C, 2)],
      redeclare Real D[size(ss1.D, 1),size(ss2.D, 2)])
      "y = G(s)*u = G(ss1)*G(ss2)*u";
  protected
    Integer nx1=size(ss1.A, 1);
    Integer nx2=size(ss2.A, 1);
algorithm
    result.A := [ss1.A,ss1.B*ss2.C; zeros(nx2, nx1),ss2.A];
    result.B := [ss1.B*ss2.D; ss2.B];
    result.C := [ss1.C,ss1.D*ss2.C];
    result.D := ss1.D*ss2.D;

end '*';

encapsulated operator function '=='
    "Check whether two linear systems have identical matrices"
    import Modelica.Math.Matrices.isEqual;
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss1 "System 1";
    input StateSpace ss2 "System 2";
    input Real eps(min=0) = 0
      "Two elements e1 and e2 of the two systems are identical if abs(e1-e2) <= eps";
    output Boolean same "=true, if the two systems are identical";
algorithm
    same := isEqual(
          ss1.A,
          ss2.A,
          eps) and isEqual(
          ss1.B,
          ss2.B,
          eps) and isEqual(
          ss1.C,
          ss2.C,
          eps) and isEqual(
          ss1.D,
          ss2.D,
          eps);
end '==';

encapsulated operator function 'String'
    "Transform state space into a String representation"
   import Modelica;
   import Modelica_LinearSystems2.StateSpace;
   import Modelica.Utilities.Strings;

   input StateSpace ss
      "State space system to be transformed in a String representation";
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

    if xNamesExist then
      Modelica.Utilities.Streams.print("xNamesExist == true");
    else
      Modelica.Utilities.Streams.print("xNamesExist == false");
    end if;

    stringMaxLength := max(size(ss.xNames, 1), min(size(ss.yNames, 1),
      11));

    if nx == 0 and sizeD == 0 then
      s := name + ".A = []\n  " + name + ".B = []\n   " + name + ".C = [] \n   " + name + ".D = []";
    else
      s := "\n" + name + ".A = \n";

  //Horizontal
  // Two alternatives when printing state names
      if xNamesExist == false then
        s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
          "x1 ";
      else
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
          ss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[
          1]), 11)) + " " + Strings.substring(
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
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
          ss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[
          1]), 11)) + " " + Strings.substring(
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
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
          ss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.xNames[
          1]), 11)) + " " + Strings.substring(
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
        s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
          ss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(ss.uNames[
          1]), 11)) + " " + Strings.substring(
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

end 'String';

encapsulated package Analysis
    "Functions to analyse state space systems represented by a StateSpace record"

  function analysis
      "Perform a system analysis based on the poles and zeros of the system"

      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Eigenvalue;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Internal;
      import Modelica.Utilities.Streams.print;

      import Modelica_LinearSystems2.Utilities.Plot;

    input StateSpace ss;

    input Internal.AnalyseOptions analyseOptions=
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
    input String fileName="eigenvalues.html"
        "Name of html-file that contains eigenvalue table";
    input String systemName="" "Name of system (used as heading in html file)";
    input String description="" "Description of system (used in html file)";
    protected
    input Boolean printStateSpaceSystem=true annotation(Dialog(enable=false));
    String dummyFileName="dummy" + fileName;
    public
    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

    protected
    input Complex j=Modelica_LinearSystems2.Math.Complex.j();
    Eigenvalue ev[size(ss.A, 1)];
    Integer nx=size(ss.A, 1);
    Integer window=0;
    Real eval[nx,2];
    Real revec[nx,nx];
    Real levec[nx,nx];
    Complex cev[size(ss.A, 1)];
    Complex systemZeros[:]=
    Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros(ss);
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
    Complex evecComplex[size(ss.A, 1),size(ss.A, 1)];
    Plot.Records.Curve curves[2];
    Plot.Records.Diagram diagram2;
    Boolean instableZeros=false;

  algorithm
    (eval,levec,revec) := Modelica_LinearSystems2.Math.Matrices.eigenValues(ss.A);

    for i in 1:nx loop
      cev[i].re := eval[i, 1];
      cev[i].im := eval[i, 2];
      ev[i].ev := cev[i];
    end for;

    (evSorted,evIndex) := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

   // Build x names
    if size(ss.xNames, 1) <> nx then
      for i in 1:nx loop
        xNames2[i] := "x[" + String(i) + "]";
      end for;
    else
      xNames2 := ss.xNames;
    end if;

  // ###### whole system checks ######
  // stability check
    isStable := true;
    for i in 1:nx loop
      isStable := isStable and ev[i].ev.re < 0;
    end for;

  //controllability check, stabilizability check
    isControllable := StateSpace.Analysis.isControllable(ss);
    isStabilizable := StateSpace.Analysis.isStabilizable(ss);

  // observability check, detectability check
    isObservable := StateSpace.Analysis.isObservable(ss);
    isDetectable := StateSpace.Analysis.isDetectable(ss);

  // analysis of single eingenvalues
    ev := StateSpace.Internal.characterizeEigenvalue(ss, ev);

  // Sort eigen values according to smallest imaginary value and restore the original order
    evSorted := Modelica_LinearSystems2.Internal.sortEigenvalue(ev);

  // analysis file

    Modelica.Utilities.Files.removeFile(fileName);
    Modelica.Utilities.Files.removeFile(dummyFileName);
    if printStateSpaceSystem then
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
      analyseOptions=analyseOptions);
    printHead1(
      ss,
      isStable,
      isControllable,
      isStabilizable,
      isObservable,
      isDetectable,
      dummyFileName,
      analyseOptions=analyseOptions);

    Modelica.Utilities.Streams.readFile(dummyFileName);

  // Plot step response
    if analyseOptions.plotStepResponse then
      Modelica.Utilities.Files.removeFile(dummyFileName);
      print("<html>\n<body>\n<p>\n<b>Step responses</b>\n</p></body></html>",
        dummyFileName);
      Modelica.Utilities.Streams.readFile(dummyFileName);
      StateSpace.Plot.step(ss=ss);
    end if;

  // Plot Bode plots
    if analyseOptions.plotFrequencyResponse then
      Modelica.Utilities.Files.removeFile(dummyFileName);
      print("<html>\n<body>\n<p>\n<b>Bode plots</b>\n</p></html>", dummyFileName);
      Modelica.Utilities.Streams.readFile(dummyFileName);
      StateSpace.Plot.bodeMIMO(ss=ss);
    end if;

   // calculate the number of real eigenvalues
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

    if analyseOptions.printEigenValues then
      if nReal > 0 then
        printHead2a(fileName, analyseOptions=analyseOptions);
        printTab1(
          evSorted,
          evIndex,
          revec,
          levec,
          nReal,
          xNames2,
          fileName,
          analyseOptions=analyseOptions);

        Modelica.Utilities.Files.removeFile(dummyFileName);
        printHead2a(dummyFileName, analyseOptions=analyseOptions);
        printTab1(
          evSorted,
          evIndex,
          revec,
          levec,
          nReal,
          xNames2,
          dummyFileName,
          analyseOptions=analyseOptions);
        Modelica.Utilities.Streams.readFile(dummyFileName);
      else
        print("<b>The system has no real eigenvalues</b><br><br>", fileName);
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print("<html><body><b>The system has no real eigenvalues</b><br><br></body></html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);

      end if;

      if nReal < nx then
        printHead2b(fileName, analyseOptions=analyseOptions);
        printTab2(
          evSorted,
          evIndex,
          revec,
          levec,
          nReal,
          xNames2,
          fileName,
          analyseOptions=analyseOptions);

        Modelica.Utilities.Files.removeFile(dummyFileName);
        printHead2b(dummyFileName, analyseOptions=analyseOptions);
        printTab2(
          evSorted,
          evIndex,
          revec,
          levec,
          nReal,
          xNames2,
          dummyFileName,
          analyseOptions=analyseOptions);
        Modelica.Utilities.Streams.readFile(dummyFileName);

    // Plot eigen values
        i := 0;
        if analyseOptions.plotEigenValues then
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
        if size(systemZeros, 1) > 0 and analyseOptions.plotInvariantZeros then
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
      else
        print("<b>The system has no conjugated complex eigenvalues</b><br><br>", fileName);
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print("<html><body><b>No conjugated complex eigenvalues</b><br><br></body></html>",
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);
      end if;

      if analyseOptions.printEigenValueProperties then
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
        printHead3(dummyFileName);
        printTab3(
          evSorted,
          evecComplex,
          evIndex,
          cev,
          nReal,
          xNames2,
          dummyFileName);
        Modelica.Utilities.Streams.readFile(dummyFileName);

      end if;
    end if;

  // ZEROS
    (zerosSorted,zerosIndex) :=
      Modelica_LinearSystems2.Math.Complex.Vectors.sortComplex(systemZeros);
    nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(zerosSorted);

    if analyseOptions.printInvariantZeros then
      if size(systemZeros, 1) > 0 then
        printHead4(fileName);
        Modelica_LinearSystems2.StateSpace.Analysis.analysis.printTab4(
          zerosSorted,
          zerosIndex,
          nReal,
          fileName);

        Modelica.Utilities.Files.removeFile(dummyFileName);
        printHead4(dummyFileName);
        printTab4(
          zerosSorted,
          zerosIndex,
          nReal,
          dummyFileName);
      else
        print("The system has no invariant zeros<br><br>", fileName);
        Modelica.Utilities.Files.removeFile(dummyFileName);
        print("<html><body><p><br><br><b>Invariant zeros</b><br>The system has no invariant zeros<br><br></body></html>",
          dummyFileName);
      end if;
      k := 0;
      for i in 1:size(systemZeros, 1) loop
        if systemZeros[i].re > 0 then
          k := k + 1;
        end if;
      end for;
      if k > 0 then
        print("<b>Note, that the system has " + String(k) + " zeros in the right complex half-plane</b>",
          fileName);
        print("<html><body><b>Note, that the system has " + String(k) + " zeros in the right complex half-plane</b></body></html>",
          dummyFileName);
      end if;

    end if;
    print("</body></html>", fileName);
    print("</body></html>", dummyFileName);
    Modelica.Utilities.Streams.readFile(dummyFileName);
    Modelica.Utilities.Files.removeFile(dummyFileName);

    print("\n\nAnalysis results have been written to file \"" +
      Modelica.Utilities.Files.fullPathName(fileName) + "\"");

  // SUB FUNCTIONS

  equation

    public
    encapsulated function printSystem
        "Print the state space system in html format on file"
        import Modelica;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;
        import Modelica_LinearSystems2;

      input StateSpace ss "state space system to analyze";
      input String fileName="systemAnalysis.html"
          "File on which the state space is written in html format";
      input String systemName="State Space System"
          "name of the state space system";
      input String description="" "Description of system (used in html file)";
      input String format=".3g" "Format of numbers (e.g. \"20.8e\")";
      protected
      Integer nx=size(ss.A, 1);
      Integer nu=size(ss.B, 2);
      Integer ny=size(ss.C, 1);
      Integer c1=integer(ceil(nx/2) - 1);
      Integer c2=integer(ceil(ny/2) - 1);
      Integer dist=8;

    algorithm
      Modelica.Utilities.Files.removeFile(fileName);
      print("<html><body><p><br><br><b>System report</b></p>", fileName);
      print("<p><br> The system <b>" + systemName + "</b><br></p>", fileName);
      print("<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
         + "cellpadding=\"3\" border=\"0\"> ", fileName);
      print("<tr><td>der(x) </td> <td>=</td> <td> Ax</td> <td> +</td><td> Bu</td></tr>
         <tr><td> y </td>     <td>=</td> <td> Cx</td> <td> + </td><td>Du</td></tr></table> <br><br>is defined by<br>",
        fileName);

    // print A and B
      print("<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
         + "cellpadding=\"3\" border=\"0\"> ", fileName);
      print("<tr><td><br></td><td><br></td><td><br></td>", fileName);
      for i1 in 1:nx loop
        print("<td style=\"text-align:center\" valign=\"top\"> " + ss.xNames[i1] + " </td>",
          fileName);
      end for;
      for i1 in 1:dist loop
        print("<td><br></td>", fileName);
      end for;
      for i1 in 1:nu loop
        print("<td>" + ss.uNames[i1] + "</td>", fileName);
      end for;

      print("</tr>", fileName);

    //print upper parts of A and B
      for i1 in 1:c1 loop
        print("<tr> <td><br></td><td><br></td><td> " + ss.xNames[i1] + " </td>",
          fileName);
        for i2 in 1:nx loop
          print("<td> " + String(ss.A[i1, i2], format=format) + " </td>", fileName);
        end for;

        for i2 in 1:dist - 1 loop
          print("<td><br></td>", fileName);
        end for;

        if nu > 0 then
          print("<td>" + ss.xNames[i1] + " </td>", fileName);
          for i2 in 1:nu loop
            print("<td> " + String(ss.B[i1, i2], format=format) + " </td>",
              fileName);
          end for;
        end if;
           //nu>0
        print("</tr>", fileName);
      end for;

    //print middle part of A and B
      print("<tr><td>A</td><td>=</td><td>" + ss.xNames[c1 + 1] + " </td>", fileName);
      for i2 in 1:nx loop
        print("<td> " + String(ss.A[c1 + 1, i2], format=format) + " </td>",
          fileName);
      end for;
      for i2 in 1:dist - 3 loop
        print("<td><br></td>", fileName);
      end for;

      if nu > 0 then
        print("<td>B</td><td>=</td><td>" + ss.xNames[c1 + 1] + " </td>", fileName);

        for i2 in 1:nu loop
          print("<td> " + String(ss.B[c1 + 1, i2], format=format) + " </td>",
            fileName);
        end for;
      else
        print("<td>B</td><td>=</td><td>[ ],  empty matrix: " +String(nx)+" by "+String(nu) + " </td>", fileName);

      end if;
           //nu>0

    //print lower parts of A and B
      for i1 in c1 + 2:nx loop
        print("<tr><td><br></td><td><br></td><td> " + ss.xNames[i1] + " </td>",
          fileName);
        for i2 in 1:nx loop
          print("<td> " + String(ss.A[i1, i2], format=format) + " </td>", fileName);
        end for;

        for i2 in 1:dist - 1 loop
          print("<td><br></td>", fileName);
        end for;
        if nu > 0 then
          print("<td>" + ss.xNames[i1] + " </td>", fileName);
          for i2 in 1:nu loop
            print("<td> " + String(ss.B[i1, i2], format=format) + " </td>",
              fileName);
          end for;
        end if;
           //nu>0
        print("</tr>", fileName);
      end for;

      print("</table>\n", fileName);

      print("<br><br>", fileName);

    // print C and D
      print("<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
         + "cellpadding=\"3\" border=\"0\">", fileName);
      print("<tr><td><br></td><td><br></td><td><br></td>", fileName);

      if ny > 0 then
        for i1 in 1:nx loop
          print("<td style=\"text-align:center\" valign=\"top\"> " + ss.xNames[i1] + " </td>",
            fileName);
        end for;
        for i1 in 1:dist loop
          print("<td><br></td>", fileName);
        end for;
        for i1 in 1:nu loop
          print("<td>" + ss.uNames[i1] + "</td>", fileName);
        end for;
      end if;
           //ny>0
      print("</tr>", fileName);

      if ny > 0 then
     //print upper parts of C and D
        for i1 in 1:c2 loop
          print("<tr> <td><br></td><td><br></td><td> " + ss.yNames[i1] + " </td>",
            fileName);
          for i2 in 1:nx loop
            print("<td> " + String(ss.C[i1, i2], format=format) + " </td>",
              fileName);
          end for;

          for i2 in 1:dist - 1 loop
            print("<td><br></td>", fileName);
          end for;
          print("<td>" + ss.yNames[i1] + " </td>", fileName);
          for i2 in 1:nu loop
            print("<td> " + String(ss.D[i1, i2], format=format) + " </td>",
              fileName);
          end for;

          print("</tr>", fileName);
        end for;
      end if;
           //ny>0

     //print middle part of C and D
      if ny > 0 then
        print("<tr><td>C</td><td>=</td><td>" + ss.yNames[c2 + 1] + " </td>",
          fileName);
        for i2 in 1:nx loop
          print("<td> " + String(ss.C[c2 + 1, i2], format=format) + " </td>",
            fileName);
        end for;
        for i2 in 1:dist - 3 loop
          print("<td><br></td>", fileName);
        end for;
        print("<td>D</td><td>=</td><td>" + (if nu>0 then ss.yNames[min(c2 + 1, ny)] else "[ ],  empty matrix: " +String(ny)+" by "+String(nu)) + " </td>",
          fileName);

        for i2 in 1:nu loop
          print("<td> " + String(ss.D[c2 + 1, i2], format=format) + " </td>",
            fileName);
        end for;
      else
        print("<tr><td>C</td><td>=</td><td>[ ],  empty matrix: " +String(ny)+" by "+String(nx) + " </td>", fileName);
        for i2 in 1:dist - 3 loop
          print("<td><br></td>", fileName);
        end for;
        print("<td>D</td><td>=</td><td>[ ],  empty matrix: " +String(ny)+" by "+String(nu) + " </td>", fileName);
      end if;
           //ny>0

     //print lower parts of C and D
      if ny > 0 then
        for i1 in c2 + 2:ny loop
          print("<tr><td><br></td><td><br></td><td> " + ss.yNames[i1] + " </td>",
            fileName);
          for i2 in 1:nx loop
            print("<td> " + String(ss.C[i1, i2], format=format) + " </td>",
              fileName);
          end for;

          for i2 in 1:dist - 1 loop
            print("<td><br></td>", fileName);
          end for;
          print("<td>" + ss.yNames[i1] + " </td>", fileName);
          for i2 in 1:nu loop
            print("<td> " + String(ss.D[i1, i2], format=format) + " </td>",
              fileName);
          end for;
          print("</tr>", fileName);
        end for;
      end if;
           // ny>0

      if description == "" then
        print("</table>\n", fileName);
      else
        print("</table>", fileName);
        print("<p>\n<b>Description</b>\n</p>", fileName);
        print(description, fileName);
      end if;

      if ny==0 and nu==0 then
        print("<p>\n<b>Note</b>, that matrices <b>B</b> and <b>C</b> are empty matrices, i.e. the system has neither inputs nor outputs!\n</p>", fileName);
      elseif ny==0 then
         print("<p>\n<b>Note</b>, that atrix <b>C</b> is empty matrix, i.e. the system has no outputs!\n</p>", fileName);
      elseif nu==0 then
         print("<p>\n<b>Note</b>, that matrix <b>B</b> is empty matrix, i.e. the system has no inputs!\n</p>", fileName);
      end if;

      print("</body></html>", fileName);

    end printSystem;

    encapsulated function printHead1
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      input Boolean isStable;
      input Boolean isControllable;
      input Boolean isStabilizable;
      input Boolean isObservable;
      input Boolean isDetectable;
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

    algorithm
      print("<html>\n<body>\n<p>\n<b>Characteristics</b>\n</p>The system\n<p>" + "</p> is ",
        fileName);
      if analyseOptions.printControllability and analyseOptions.printObservability then
        print((if isStable then " " else "not ") + "stable" + "\n<br>" + (if 
          isStable then if isControllable then "and it is " else "but it is not " else 
                if isControllable then "but it is " else "and it is not ") + "controllable"
           + (if isStable then "" else "\n<br>" + (if isControllable then " and therefore it is " else 
                if isStabilizable then " but it is " else "and is not ") + "stabilizable")
           + "\n<br> The system is " + (if isObservable then " " else "not ") + "observable"
           + (if isStable then "" else "\n<br>" + (if isObservable then " and therefore it is " else 
                if isDetectable then " but it is " else "and is not ") + "detectable")
           + "\n<br></br>", fileName);
      elseif not analyseOptions.printObservability and analyseOptions.printControllability then
        print((if isStable then " " else "not ") + "stable" + "\n<br>" + (if 
          isStable then if isControllable then "and it is " else "but it is not " else 
                if isControllable then "but it is " else "and it is not ") + "controllable"
           + (if isStable then "" else "\n<br>" + (if isControllable then " and therefore it is " else 
                if isStabilizable then " but it is " else "and is not ") + "stabilizable")
           + "\n<br></br>", fileName);
      elseif not analyseOptions.printControllability and analyseOptions.printObservability then
        print((if isStable then " " else "not ") + "stable" + "\n<br> The system is "
           + (if isObservable then " " else "not ") + "observable" + (if isStable then 
                "" else "\n<br>" + (if isObservable then " and therefore it is " else 
                if isDetectable then " but it is " else "and is not ") + "detectable")
           + "\n<br></br>", fileName);
      else
        print((if isStable then " " else "not ") + "stable" + "\n<br></br>",
          fileName);
      end if;
      print("</body></html>", fileName);

    end printHead1;

    public
    encapsulated function printHead2a
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

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

    algorithm
      if analyseOptions.printEigenValueProperties then
        print("<html>\n<body>\n<b><big>Eigenvalues analysis</big></b><br><br><b>Real eigenvalues</b>\n<br>"
           + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
           + "cellpadding=\"3\" border=\"1\">\n" + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           + "<td> number </td><td> eigenvalue </td> <td> T [s] </td>  <td> characteristics </td><td> contribution to states</td></tr>",
          fileName);
      else
        print("<html>\n<body>\n<b><big>Eigenvalues analysis</big></b><br><br><b>Real eigenvalues</b>\n<br>"
           + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
           + "cellpadding=\"3\" border=\"1\">\n" + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           + "<td> number </td><td> eigenvalue </td> <td> T [s] </td>  <td> characteristics </td></tr>",
          fileName);
      end if;
    end printHead2a;

    public
    encapsulated function printHead2b
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

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

    algorithm
      print("<html>\n<body>", fileName);
      if analyseOptions.printEigenValueProperties then
        print("<b>Conjugated complex pairs of eigenvalues</b>\n<br>" + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
           + "cellpadding=\"3\" border=\"1\">\n" + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           + "<td> number </td> <td> eigenvalue </td><td> freq. [Hz] </td> <td> damping </td><td> characteristics </td>  <td> contribution to states</td></tr>",
          fileName);
      else
        print("<b>Conjugated complex pairs of eigenvalues</b>\n<br>" + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
           + "cellpadding=\"3\" border=\"1\">\n" + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
           + "<td> number </td> <td> eigenvalue </td><td> freq. [Hz] </td> <td> damping </td><td> characteristics </td> </tr>",
          fileName);
      end if;
    end printHead2b;

    encapsulated function printHead3
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

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

    algorithm
      print("<html><body>", fileName);
      print("<p> In the tables above, the column <b>contribution to states</b> lists for each eigenvalue the states to which the"
         + "corresponding modal state contributes most. This information is based on the "
         + "two largest absolute values of the corresponding right eigenvector (if the second large value "
         + "is less than 5 % of the largest contribution, it is not shown).<br>" + "</p><p>"
         + "In the next table, for each state in the column <b>correlation to modal states</b>, the modal  "
         + "states which contribute most to the coresponding state are summarized, i.e. the state is mostly composed of these modal states "
         + "This information is based on the two largest absolute values of row i of the "
         + "eigenvector matrix that is associated with eigenvalue i (if the second large value  "
         + "is less than 5 % of the largest contribution, it is not shown). This only holds "
         + "if the modal states are in the same order of magnitude. Otherwise, the modal states "
         + "listed in the last column might be not the most relevant one. </p><br><br> "
         + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
         + "cellpadding=\"3\" border=\"1\"> " + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
         + "<td> state </td> <td> composition </td> <td> eigenvalue #</td> <td> freq. [Hz] </td> <td> damping </td>  </td> <td> T [s] </td></tr>",
        fileName);

    end printHead3;

    encapsulated function printHead4
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;

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

    algorithm
      print("<html><body>", fileName);
      print("<p><br><br><b>Invariant zeros</b><br>" + "<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" "
         + "cellpadding=\"3\" border=\"1\">" + "<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">"
         + "<td> number </td> <td> invariant zero </td><td> Time constant [s] </td> <td> freq. [Hz] </td> <td> damping </td></tr>",
        fileName);
    end printHead4;

    encapsulated function printTab1
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Modelica_LinearSystems2.Math.Complex;

      input Eigenvalue evSorted[:];
      input Integer evIndex[size(evSorted, 1)];
      input Real r_evec[size(evSorted, 1),size(evSorted, 1)];
      input Real l_evec[size(evSorted, 1),size(evSorted, 1)];
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
        if analyseOptions.printEigenValueProperties then
          print("<tr>\n <td style=\"text-align:center\"> " + number + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].timeConstant, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + (if evSorted[i].isStable then "" else "not ") + "stable, " + (if 
            evSorted[i].isStable then (if evSorted[i].isControllable then "" else 
                  "not ") + "controllable, " else (if evSorted[i].isStabilizable then 
                  "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable then 
                  (if evSorted[i].isObservable then "" else "not ") + "observable " else 
                  (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td> <td style=\"text-align:left\"> &nbsp; " + " z[" + String(i)
             + "]" + " contributes to " + xNames2[r_maxIndex1] + " with " +
            String(r_absMax1, format=".3g") + " %<br>" + (if r_two then "&nbsp; " + " z["
             + String(i) + "]" + " contributes to " + xNames2[r_maxIndex2] + " with "
             + String(r_absMax2, format=".3g") + " %" else "") + " </td> </tr> ",
            fileName);

        else
          print("<tr>\n <td style=\"text-align:center\"> " + number + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].timeConstant, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + (if evSorted[i].isStable then "" else "not ") + "stable, " + (if 
            evSorted[i].isStable then (if evSorted[i].isControllable then "" else 
                  "not ") + "controllable, " else (if evSorted[i].isStabilizable then 
                  "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable then 
                  (if evSorted[i].isObservable then "" else "not ") + "observable " else 
                  (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td> </tr> ", fileName);
        end if;
        i := j;
      end while;

      print("</table><br><br>\n\n</body></html>", fileName);
    end printTab1;

    encapsulated function printTab2
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Modelica_LinearSystems2.Math.Complex;

      input Eigenvalue evSorted[:];
      input Integer evIndex[size(evSorted, 1)];
      input Real r_evec[size(evSorted, 1),size(evSorted, 1)];
      input Real l_evec[size(evSorted, 1),size(evSorted, 1)];
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
        if analyseOptions.printEigenValueProperties then
          print("<tr>\n <td style=\"text-align:left\"> " + number + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " &plusmn; " + String(
            evSorted[i].ev.im, format="12.4e") + "j" + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].frequency, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].damping, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + (if evSorted[i].isStable then "" else "not ") + "stable, " + (if 
            evSorted[i].isStable then (if evSorted[i].isControllable then "" else 
                  "not ") + "controllable, " else (if evSorted[i].isStabilizable then 
                  "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable then 
                  (if evSorted[i].isObservable then "" else "not ") + "observable " else 
                  (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td> <td style=\"text-align:left\"> &nbsp; " + " z[" + number + "]"
             + " contribute to " + xNames2[r_maxIndex1] + " with " + String(
            r_absMax1, format=".3g") + " %<br>" + (if r_two then "&nbsp; " + " z["
             + number + "]" + " contribute to " + xNames2[r_maxIndex2] + " with " +
            String(r_absMax2, format=".3g") + " %" else "") + " </td> </tr> ",
            fileName);
        else
          print("<tr>\n <td style=\"text-align:left\"> " + number + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].ev.re, format="14.4e") + " &plusmn; " + String(
            evSorted[i].ev.im, format="12.4e") + "j" + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].frequency, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + String(evSorted[i].damping, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; "
             + (if evSorted[i].isStable then "" else "not ") + "stable, " + (if 
            evSorted[i].isStable then (if evSorted[i].isControllable then "" else 
                  "not ") + "controllable, " else (if evSorted[i].isStabilizable then 
                  "" else "not ") + "stabilizable, ") + (if evSorted[i].isStable then 
                  (if evSorted[i].isObservable then "" else "not ") + "observable " else 
                  (if evSorted[i].isDetectable then "" else "not ") + "detectable ")
             + " </td> </tr> ", fileName);
        end if;
        i := j;
      end while;

      print("</table>", fileName);
      print("</table><br><br>\n\n</body></html>", fileName);
    end printTab2;

    encapsulated function printTab3
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Modelica_LinearSystems2.Math.Complex;

      input Eigenvalue evSorted[:];
      input Complex evecComplex[:,:];
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
        v_normalized := Complex.Vectors.normalize(evecComplex[i, :]);
        first := true;
        two := false;
        absMax1 := 0;
        maxIndex1 := 0;
        absMax2 := 0;
        maxIndex2 := 0;
        j := 1;
        abs_v_normalized := Complex.Vectors.norm(v_normalized, 1);
        while j <= nx loop
          if cev[j].im == 0 then
            v := abs(v_normalized[j].re);
            k := j;
            j := j + 1;
          else
            v := 2*Complex.'abs'(v_normalized[j]);
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
        absMax1 := absMax1/abs_v_normalized;
        absMax2 := absMax2/abs_v_normalized;

        if absMax2 < 0.05*absMax1 then
          two := false;
        end if;

         // Determine frequency and number of corresponding eigenvalue
        (w1,d1) := Complex.frequency(cev[maxIndex1]);
        iw1 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex1, evIndex);
        if iw1 <= nReal then
          number1 := String(iw1);
        else
          number1 := String(iw1) + "/" + String(iw1 + 1);
        end if;

        if two then
          (w2,d2) := Complex.frequency(cev[maxIndex2]);
          iw2 := Modelica_LinearSystems2.Math.Vectors.find(maxIndex2, evIndex);
          if iw2 <= nReal then
            number2 := String(iw2);
          else
            number2 := String(iw2) + "/" + String(iw2 + 1);
          end if;
        end if;

        print("<tr>\n <td style=\"text-align:left\"> &nbsp; " + xNames2[i] + " </td> <td style=\"text-align:left\"> &nbsp; "
           + " is composed of " + String(100*absMax1, format="5.1f") + "% by z[" +
          number1 + "]" + (if two then " <br>" + " &nbsp; " + " is composed of " +
          String(100*absMax2, format="5.1f") + "% by z[" + number2 + "]" else "") + " </td> <td style=\"text-align:center\"> &nbsp; "
           + number1 + (if two then "<br> &nbsp; " + number2 else Strings.repeat(9))
           + " </td> <td style=\"text-align:center\"> &nbsp; " + (if iw1 <= nReal then 
                "---" else String(w1, format="9.4f")) + (if two then "<br> &nbsp; "
           + (if iw2 <= nReal then "---" else String(w2, format="9.4f")) else 
          Strings.repeat(9)) + " </td> <td style=\"text-align:center\"> &nbsp; " +
          (if iw1 <= nReal then "---" else String(d1, format="9.4f")) + (if two then 
                "<br> &nbsp; " + (if iw2 <= nReal then "---" else String(d2,
          format="9.4f")) else "") + " </td> <td style=\"text-align:center\"> &nbsp; "
           + (if (iw1 <= nReal) then String(evSorted[i].timeConstant, format="9.4f") else 
                "---") + (if two then "<br> &nbsp; " + (if (iw2 <= nReal and abs(
          cev[maxIndex2].re) > 1e-10) then String(1/abs(cev[maxIndex2].re),
          format="9.4f") else "---") else "") + " </td> </tr> ", fileName);

      end for;
      print("</table></body></html>", fileName);

    end printTab3;

    encapsulated function printTab4
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica_LinearSystems2;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Modelica_LinearSystems2.Math.Complex;

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
        timeConstant := if abs(systemZeros[i].re) > 10*Modelica.Constants.eps then 
                1/abs(systemZeros[i].re) else 1/(10*Modelica.Constants.eps);

        print("<tr>\n <td style=\"text-align:left\"> &nbsp; " + number + " </td> <td> &nbsp; "
           + String(systemZeros[i].re, format="14.4e") + " </td> <td> &nbsp; " +
          String(timeConstant, format="9.4f") + " </td> <td style=\"text-align:center\"> &nbsp; "
           + "---" + " </td> <td style=\"text-align:center\"> &nbsp; " + "---" + " </td> </tr> ",
          fileName);

      end for;

      for i in nReal + 1:2:nz loop
        number := String(i) + "/" + String(i + 1);
        number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;

       // Determine frequency and number of corresponding zero
        (freq,damp) := Complex.frequency(systemZeros[i]);

        print("<tr>\n <td style=\"text-align:left\"> &nbsp; " + number + " </td> <td style=\"text-align:left\"> &nbsp; "
           + String(systemZeros[i].re, format="14.4e") + " &plusmn; " + String(
          systemZeros[i].im, format="12.4e") + "j" + " </td> <td style=\"text-align:center\"> &nbsp; "
           + "---" + " </td> <td style=\"text-align:left\"> &nbsp; " + String(
          freq, format="9.4f") + " </td> <td style=\"text-align:left\"> &nbsp; " +
          String(damp, format="9.4f") + " </td> </tr> ", fileName);

      end for;

      print("</table><br><br><br>", fileName);
    end printTab4;
    annotation (interactive=true, Documentation(info="<html>
<h4><font style=\"color: #008000; \">Syntax</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(ss);
   or
Modelica_LinearSystems2.StateSpace.Analysis.<b>analysis</b>(ss, analyseOptions=<a href=\"Modelica://Modelica_LinearSystems2.Internal.AnalyseOptions\">analyseOptions</a>, fileName, systemName, description);
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>

This function analyzes a state space system <br>
<pre>    der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
                                <label for=\"eqn1\">(1)</label>
        <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>
        <b>x</b>(t=0) = <b>x</b><sub>0</sub>

</pre>
based on its poles, i.e. the eigenvalues, and the zeros of the system.
The system will be checked for stability, controllability and observability. In the case that the system is not stable stabilizability and detectability are examined. Furthermore, stability, controllability, observability, stabilizability, and detectability are indicated for each eigenvalue.
<br>

<br><br>
<b>Stability</b><br>

System (1) is stable if and only if all eigenvalues of the matrix <b>A</b> have negative real parts.
<br>
The calculation of the eigenvalues is based on the LAPACK routine dgeev.
<br><br>
<b>Controllability</b><br>
System (1) is said to be controllable if, starting from any initial state <b>x</b><sub>0</sub>, the system can be driven by appropriate inputs to any final state <b>x</b><sub>1</sub> within some finite time window. Equivalent is that the eigenvalues of <b>A</b>-<b>BK</b> can  arbitrarily be assigned by an appropriate choice of the matrix <b>K</b>.
<br>
<br>
<b>Stabilizability</b><br>

System (1) is said to be stabilizable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are controllable. Therefore, a controllable system is always stabilizable. An equivalent definition of stabilizability is, that a system is said to be stabilizable if there exist a matrix <b>K</b> such that <b>A</b>-<b>BK</b> is stable.
<br>
<br>
<b>Observability</b><br>

System (1) is said to be observable if the (arbitrary) initial state <b>x</b><sub>0</sub> can be uniquely determined from any state <b>x</b>(t<sub>1</sub>), t<sub>1</sub>>0, from the knowledge of the input <b>u</b>(t) and output <b>y</b>(t). With other words,  from the system's outputs it is possible to determine the behavior of the entire system. Equivalent is, that the eigenvalues of <b>A</b>-<b>LC</b> can be arbitrarily be assigned by an appropriate choice of matrix <b>L</b>.<br>
Observability is called the dual concept of controllability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is observable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is controllable.

<br>
<br>
<b>Detectability</b><br>

System (1) is said to be detectable if all the unstable eigenvalues, i.e. all <tt>s</tt> with Re(<tt>s</tt>)>=0, of <b>A</b> are observable. Therefore, a observable system is always detectable. An equivalent definition of detectability is, that a system is said to be detectable if there exist a matrix <b>L</b> such that <b>A</b>-<b>LC</b> is stable.
Detectability is called the dual concept of stabilizability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is detectable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is stabilizable.<br>
<br>
<b>Algorithm to test controllability/stabilizability and observability/detectability respectively</b> <br>

The test of controllability and stabilizability is performed with the staircase algorithm which transforms the system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) into the controller-Hessenberg form (<b>A</b><sub>H</sub>, <b>B</b><sub>H</sub>, <b>C</b><sub>H</sub>, <b>D</b>) with <b>A</b><sub>H</sub> is a block upper Hessenberg matrix and <b>B</b><sub>H</sub>=[<b>B</b><sub>1</sub>; 0] with triangular matrix <b>B</b><sub>1</sub> with rank(<b>B</b><sub>1</sub>) = rank(<b>B</b>).<br>
In <b>A</b><sub>H</sub>=[<b>A</b><sub>c</sub>, *,0, <b>A</b><sub>nc</sub>) the eigenvalues of the matrices <b>A</b><sub>c</sub> and <b>A</b><sub>nc</sub> are the controllable eigenvalues and uncontrollable eigenvalues of <b>A</b> respectively.<br>
The test of observability and detectability is performed by testing the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) with respect to controllability and stabilizability.<br>
<br>
<b>Solution of a linear time invariant system </b><br>

The solution <b>x</b>(t) of the initial value problem (1) consists of the homogeneous part (zero input response) <b>x</b><sub>h</sub>(t) and the inhomogeneous part x<sub>i</sub>(t). The zero input solution is given by
<pre>
<b>x</b><sub>h</sub>(t) = exp(<b>A</b>*(t-t<sub>0</sub>))<b>x</b><sub>0</sub>.
</pre>
The system can also be represented as a linear combination of the modal states <b>z</b>, 
<pre>
<b>x</b> = <b>V</b><b>z</b>
</pre>
i.e. the states of a similar system, with
<pre>
der(<b>z</b>) = <b>V</b><sup>-1</sup><b>AVz</b> + <b>V</b><sup>-1</sup><b>B</b><b>u</b>
</pre>
where the system matrix <b>V</b><sup>-1</sup><b>AV</b> is the real Jordan form. For single real eigenvectors the system is decoupled, i.e. the solution of the modal states are denoted by
<pre>
z<sub>i</sub> = exp(s<sub>i</sub> t)*z<sub>0i</sub>

</pre>

The behavior of the modal states is determined as the solution of a linear first order differential equation for real eigenvalues. Since this behavior is well known, the behavior of the x<sub>i</sub> can at least roughly be estimated by means of the behavior of the most relevant modal states. Therefore, the contribution of the modal states to the states is computed as an indication of the original system behavior.
<br>
<p>
<b>Contribution of the modal states to the states</b><br>
Generally, as described above, the states of the system can be described as linear combination of modal states and, therefore, the states can be characterized to a certain extend by the modal states if the proportions of the combination are known. Hence, for each modal state z<sub>i</sub> of the vector <b>z</b> the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub>| of the corresponding right eigenvector <b>v</b><sub>i</sub> indicate the proportion of <b>z</b><sub>i</sub> that is contributed to the state x<sub>j</sub>.<br>
On the other hand, the composition of xi is indicated by the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>|, i.e. the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>| of the corresponding row <b>v</b><sub>i</sub><sup>T</sup> of the eigenvector matrix <b>V</b> indicate the proportion of the state x<sub>i</sub> that is contributed by the modal state z<sub>j</sub>.

</p>



<h4><font color=\"#008000\">Example</font></h4>
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



<body>
<p>
<b>System report</b>
</p>
<p><br> The system <b>Demonstation System</b>
</p>
<br><table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"0\">
 
<tr><td>der(x) </td><td>=</td><td> Ax</td><td> +</td><td> Bu</tr><td> y </td><td>=</td><td> Cx</td><td> + Du</td></tr></table> <br>is defined by<br>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"0\">
 
<tr><td><br></td><td><br></td><td><br></td>
<td style=\"text-align:center\" valign=\"top\"> x1 </td>
<td style=\"text-align:center\" valign=\"top\"> x2 </td>
<td style=\"text-align:center\" valign=\"top\"> x3 </td>
<td style=\"text-align:center\" valign=\"top\"> x4 </td>
<td style=\"text-align:center\" valign=\"top\"> x5 </td>
<td style=\"text-align:center\" valign=\"top\"> x6 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>u1</td>
<td>u2</td>
</tr>
<tr> <td><br></td><td><br></td><td> x1 </td>
<td> -3 </td>
<td> 2 </td>
<td> -3 </td>
<td> 4 </td>
<td> 5 </td>
<td> 6 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>x1 </td>
<td> 1 </td>
<td> 0 </td>
</tr>
<tr> <td><br></td><td><br></td><td> x2 </td>
<td> 0 </td>
<td> 6 </td>
<td> 7 </td>
<td> 8 </td>
<td> 9 </td>
<td> 4 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>x2 </td>
<td> 0 </td>
<td> 1 </td>
</tr>
<tr><td>A</td><td>=</td><td>x3 </td>
<td> 0 </td>
<td> 2 </td>
<td> 3 </td>
<td> 0 </td>
<td> 78 </td>
<td> 6 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>B</td><td>=</td><td>x3 </td>
<td> 1 </td>
<td> 0 </td>
<tr><td><br></td><td><br></td><td> x4 </td>
<td> 0 </td>
<td> 1 </td>
<td> 2 </td>
<td> 2 </td>
<td> 3 </td>
<td> 3 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>x4 </td>
<td> 0 </td>
<td> 1 </td>
</tr>
<tr><td><br></td><td><br></td><td> x5 </td>
<td> 0 </td>
<td> 13 </td>
<td> 34 </td>
<td> 0 </td>
<td> 0 </td>
<td> 1 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>x5 </td>
<td> 1 </td>
<td> 0 </td>
</tr>
<tr><td><br></td><td><br></td><td> x6 </td>
<td> 0 </td>
<td> 0 </td>
<td> 0 </td>
<td> -17 </td>
<td> 0 </td>
<td> 0 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>x6 </td>
<td> 0 </td>
<td> 1 </td>
</tr>
</table>
<p>
<br><br>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"0\">
 
<tr><td><br></td><td><br></td><td><br></td>
<td style=\"text-align:center\" valign=\"top\"> x1 </td>
<td style=\"text-align:center\" valign=\"top\"> x2 </td>
<td style=\"text-align:center\" valign=\"top\"> x3 </td>
<td style=\"text-align:center\" valign=\"top\"> x4 </td>
<td style=\"text-align:center\" valign=\"top\"> x5 </td>
<td style=\"text-align:center\" valign=\"top\"> x6 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>u1</td>
<td>u2</td>
</tr>
<tr><td>C</td><td>=</td><td>y1 </td>
<td> 0 </td>
<td> 0 </td>
<td> 1 </td>
<td> 0 </td>
<td> 1 </td>
<td> 0 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>D</td><td>=</td><td>y1 </td>
<td> 0 </td>
<td> 0 </td>
<tr><td><br></td><td><br></td><td> y2 </td>
<td> 0 </td>
<td> 1 </td>
<td> 0 </td>
<td> 0 </td>
<td> 1 </td>
<td> 1 </td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td><br></td>
<td>y2 </td>
<td> 0 </td>
<td> 0 </td>
</tr>
</table>
<p>
</body>
<body>
<p>
<b>Description</b>
</p>
System to demonstrate the usage of Modelica_LinearSystems2.StateSpace.Analysis.anlysis()
</body>
<body>
<p>
<b>Characteristics</b>
</p>The system
<p></p> is 
not stable
<br>but it is controllable
<br> and therefore it is stabilizable
<br> The system is not observable
<br> but it is detectable
<br></br>
<b><big>Eigenvalues analysis</big></b><br><br><b>Real eigenvalues</b>
<br><table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"1\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td><td> eigenvalue </td> <td> T [s] </td>  <td> characteristics </td><td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:center\">       1 </td> <td style=\"text-align:left\"> &nbsp;   -4.9874e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0201 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, observable  </td> <td style=\"text-align:left\"> &nbsp;  z[1] contributes to x3 with 54.6 %<br>&nbsp;  z[1] contributes to x5 with 37 % </td> </tr> 
<tr>
 <td style=\"text-align:center\">       2 </td> <td style=\"text-align:left\"> &nbsp;   -3.0000e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3333 </td> <td style=\"text-align:left\"> &nbsp; stable, controllable, not observable  </td> <td style=\"text-align:left\"> &nbsp;  z[2] contributes to x1 with 100 %<br> </td> </tr> 
<tr>
 <td style=\"text-align:center\">       3 </td> <td style=\"text-align:left\"> &nbsp;    2.9891e+000 </td> <td style=\"text-align:left\"> &nbsp;    0.3346 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[3] contributes to x2 with 51.9 %<br>&nbsp;  z[3] contributes to x1 with 23.9 % </td> </tr> 
<tr>
 <td style=\"text-align:center\">       4 </td> <td style=\"text-align:left\"> &nbsp;    5.5825e+001 </td> <td style=\"text-align:left\"> &nbsp;    0.0179 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[4] contributes to x3 with 48.4 %<br>&nbsp;  z[4] contributes to x5 with 32.5 % </td> </tr> 
</table><br><br>


<b>Conjugated complex pairs of eigenvalues</b>
<br><table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"1\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> eigenvalue </td><td> freq. [Hz] </td> <td> damping </td><td> characteristics </td>  <td> contribution to states</td></tr>
<tr>
 <td style=\"text-align:left\">     5/6 </td> <td style=\"text-align:left\"> &nbsp;    1.0299e+000 &plusmn;  6.5528e+000j </td> <td style=\"text-align:left\"> &nbsp;    1.0557 </td> <td style=\"text-align:left\">> &nbsp;   -0.1553 </td> <td style=\"text-align:left\"> &nbsp; not stable, stabilizable, detectable  </td> <td style=\"text-align:left\"> &nbsp;  z[    5/6] contribute to x6 with 35.9 %<br>&nbsp;  z[    5/6] contribute to x2 with 20.6 % </td> </tr> 
</table>
<p>
In the tabel above, the column <b>contribution to states</b> lists for each eigenvalue the states to which thecorresponding modal state contributes most. This information is based on the
two largest absolute values of the corresponding right eigenvector (if the second large value 
is less than 5 % of the largest contribution, it is not shown).
<br> 
</p>
<p>
In the next table, for each state in the column <b>correlation to modal states</b>, the modal 
states which contribute most to the coresponding state are summarized, i.e. the state is mostly composed of these modal states
This information is based on the two largest absolute values of row i of the
eigenvector matrix that is associated with eigenvalue i (if the second large value 
is less than 5 % of the largest contribution, it is not shown). This only holds
if the modal states are in the same order of magnitude. Otherwise, the modal states
listed in the last column might be not the most relevant one.
</p>
<table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"1\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> state </td> <td> composition </td> <td> eigenvalue #</td> <td> freq. [Hz] </td> <td> damping </td>  </td> <td> T [s] </td></tr>
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
<br><br><b>Invariant zeros</b><br><table style=\"font-size:10pt; font-family:Arial; border-collapse:collapse; text-align:right\" cellpadding=\"3\" border=\"1\">
<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\"><td> number </td> <td> invariant zero </td><td> Time constant [s] </td> <td> freq. [Hz] </td> <td> damping </td></tr>
<tr>
 <td style=\"text-align:left\"> &nbsp;       1 </td> <td> &nbsp;   -5.4983e+001 </td> <td> &nbsp;    0.0182 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr> 
<tr>
 <td style=\"text-align:left\"> &nbsp;       2 </td> <td> &nbsp;   -3.0000e+000 </td> <td> &nbsp;    0.3333 </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:center\"> &nbsp; --- </td> </tr> 
<tr>
 <td style=\"text-align:left\"> &nbsp;     3/4 </td> <td style=\"text-align:left\"> &nbsp;    3.2417e+000 &plusmn;  5.6548e+000j </td> <td style=\"text-align:center\"> &nbsp; --- </td> <td style=\"text-align:left\"> &nbsp;    1.0374 </td> <td style=\"text-align:left\"> &nbsp;   -0.4973 </td> </tr> 
</table>
</body>
</html>
"),             Documentation(info="<html>
 
Function <b>Modelica_LinearSystems2.StateSpace.Analysis.analysis</b> analyzes a state space system <br>
<pre>    der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
                                <label for=\"eqn1\">(1)</label>
        <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>  
        <b>x</b>(t=0) = <b>x</b><sub>0</sub> 
</pre>
based on its poles, i.e. the eigenvalues, and the zeros of the system.
The system will be checked for stability, controllability and observability. If the case that the system is not stable stabilizability and detectability are examined. Furthermore, stability, controllability, observability, stabilizability, and detectability are indicated for each pole.
<br>
 
Generally, The states of the system can be described as linear combination of modal states (see below) and, therefore, the states can be characterized to a certain extend by the
modal states if the proportions of the combination are known. Hence, for each modal state z<sub>i</sub> of the vector <b>z</b> the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub>|$ of the corresponding right eigenvector <b>v</b><sub>i</sub> indicates the proportion of <b>z</b><sub>i</sub> that is contributed to the state x<sub>j</sub>.<br>
On the other hand, the composition of xi is indicated by the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>|, i.e. the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>| of the corresponding row <b>v</b><sub>i</sub><sup>T</sup> of the eigenvector matrix <b>V</b> indicates the proportion of the state x<sub>i</sub> that is contributed by the modal state z<sub>j</sub>. 
 
<br><br>
<b>Stability</b><br>
 
System (1) is stable if and only if all eigenvalues of the matrix <b>A</b> have negative real parts.
<br>
The calculation of the eigenvalues and also of the eigenvalues uses LAPACK routine dgeev.
<br><br>
<b>Controllability</b><br>
System (1) is said to be controllable if starting from any initial state <b>x</b>0, the system can be driven by appropriate inputs to any final state <b>x</b>1 within some finite time window. Equivalent is, that the eigenvalues of <b>A</b>-<b>BK</b> can be arbitrarily be assigned by an appropriate choice of matrix <b>K</b>. 
<br>
<br>
<b>Stabilizability</b><br>
 
System (1) is said to be stabilizable if all the unstable eigenvalues, i.e. all &lambda with Re(&lambda)>=0, of <b>A</b> are controllable. Therefore, a controllable system is always stabilizable. An equivalent definition of stabilizability is, that a system is said to be stabilizable if there exist a matrix K such that A-BK is stable.
<br>
<br>
<b>Observability</b><br>
 
System (1) is said to be observable if the (arbitrary) initial state x0 can be uniquely determined from any sate x(t1), t1>0, from the knowledge of the input u(t) and output y(t). With other words,  from the system's outputs it is possible to determine the behaviour of the entire system. Equivalent is, that the eigenvalues of <b>A</b>-<b>LC</b> can be arbitrarily be assigned by an appropriate choice of matrix <b>L</b>.<br>
Observability is called the dual concept of controllability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is observable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is controllable.
 
<br>
<br>
<b>Detectability</b><br>
 
System (1) is said to be detectable if all the unstable eigenvalues, i.e. all &lambda with Re(&lambda)>=0, of <b>A</b> are observable. Therefore, a observable system is always detectable. An equivalent definition of detectability is, that a system is said to be detectable if there exist a matrix <b>L</b> such that <b>A</b>-<b>LC</b> is stable.
Detectability is called the dual concept of stabilizability, since a system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) is detectable if the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) is stabilizable.<br>
<br>
<b>Algorithm to test controllability/stabilizability and observability/detectability respectively</b> <br>
 
The test of controllability and stabilizability is performed with the staircase algorithm which transforms the system (<b>A</b>,<b>B</b>,<b>C</b>,<b>D</b>) into the controller-Hessenberg form (<b>A</b><sub>H</sub>, <b>B</b><sub>H</sub>, <b>C</b>H, <b>D</b>) with <b>A</b>H is a block upper Hessenberg matrix and <b>B</b><sub>H</sub>=[<b>B</b>1; 0] with triangular matrix <b>B</b>1 with rank(<b>B</b>1) = rank(<b>B</b>).<br>
In <b>A</b><sub>H</sub>=[<b>A</b>c, *,0, <b>A</b>nc) the eigenvalues of the matrices <b>A</b>c and <b>A</b>nc are the controllable eigenvalues and uncontrollable eigenvalues of <b>A</b> respectively.<br>
The test of observability and detectability is performed by testing the system (<b>A</b><sup>T</sup>, <b>C</b><sup>T</sup>, <b>B</b><sup>T</sup>, <b>D</b><sup>T</sup>) with respect to controllability and stabilizability.<br>
<br>
<b>Solution of a linear time invariant system </b><br>
 
The solution x(t) of the initial value problem (1) consists of the homogeneous part (zero input response) x<sub>h</sub>(t) and the inhomogeneous part x<sub>i</sub>(t). The zero input solution is given by
<pre>
<b>x</b><sub>h</sub>(t) = exp(<b>A</b>*(t-t<sub>0</sub>))<b>x</b><sub>0</sub>.
</pre>
The system can also be represented as a linear combination of the modal system, i.e. the solution a similar system
<pre>
<b>x</b> = <b>V</b><b>z</b>
</pre>
with
<pre>
der(z) = V-1AVz + V-1Bu
</pre>
and the real Jordan form  V-1AV. For single real eigenvectors is decoupled, i.e. the solution of the modal states are denoted ba
<pre>
z_i = exp(&lambda_i t)*z0i
</pre>
 
The behavior of the modal states is determined as the solution of a linear first order differential equation for real eigenvalues. Since this behavior is well known, the behavior of the xi can at least roughly be estimated by means of the behavior of the most relevant modal states. Therefore, the contribution of the modal states to the states is computed an .
<br>
 
<b>Contribution of the modal states to the states</b><br>
Generally, as described above, the states of the system can be described as linear combination of modal states and, therefore, the states can be characterized to a certain extend by the
modal states if the proportions of the combination are known. Hence, for each modal state z<sub>i</sub> of the vector <b>z</b> the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub>|$ of the corresponding right eigenvector <b>v</b><sub>i</sub> indicates the proportion of <b>z</b><sub>i</sub> that is contributed to the state x<sub>j</sub>.<br>
On the other hand, the composition of xi is indicated by the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>|, i.e. the elements |v<sub>i,j</sub>|/|<b>v</b><sub>i</sub><sup>T</sup>| of the corresponding row <b>v</b><sub>i</sub><sup>T</sup> of the eigenvector matrix <b>V</b> indicates the proportion of the state x<sub>i</sub> that is contributed by the modal state z<sub>j</sub>. 
 
 
 
<a href=\"file:ls2Docu2.pdf\">text.pdf</a>
 
 
 
 
</html>"));
  end analysis;

 encapsulated function timeResponse
      "Calculate the time response of a state space system"

   import Modelica;
   import Modelica_LinearSystems2;
   import Modelica_LinearSystems2.StateSpace;
   import Modelica_LinearSystems2.Types.TimeResponse;

   input TimeResponse response=TimeResponse.Step;
   extends Modelica_LinearSystems2.Internal.timeResponseMask2(redeclare Real y[:,size(sc.C, 1),if response == TimeResponse.Initial then 1 else size(sc.B, 2)],
       redeclare Real x_continuous[:,size(sc.A, 1),if response == TimeResponse.Initial then 1 else size(sc.B, 2)]);  // Input/Output declarations of time response functions

   input Real x0[size(sc.A, 1)]=zeros(size(sc.A, 1)) "Initial state vector";

    protected
   Real dtVar;
   Real tSpanVar;
   Integer samples;
   Real u[:,size(sc.B, 2)];
   Real new_x[size(sc.A, 1),1];
   Real x[size(sc.A, 1),1]=zeros(size(sc.A, 1), 1);
   Modelica_LinearSystems2.DiscreteStateSpace sd(
     redeclare Real A[size(sc.A, 1),size(sc.A, 2)],
     redeclare Real B[size(sc.B, 1),size(sc.B, 2)],
     redeclare Real C[size(sc.C, 1),size(sc.C, 2)],
     redeclare Real D[size(sc.D, 1),size(sc.D, 2)],
     redeclare Real B2[size(sc.B, 1),size(sc.B, 2)]);
   Real i1;
   Real i2;

 algorithm
      // set sample time and simulation time span
  if (dt == 0 and tSpan == 0) then
    (dtVar,tSpanVar) := Modelica_LinearSystems2.Internal.timeResponseSamples(
      sc);
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
  y := if response == TimeResponse.Initial then zeros(samples, size(sc.C, 1),1) else zeros(samples, size(sc.C, 1),size(sc.B, 2));
  x_continuous :=  if response == TimeResponse.Initial then zeros(samples, size(sc.A, 1),  1) else zeros(samples, size(sc.A, 1),  size(sc.B, 2));

 if response == TimeResponse.Initial then
      sd := Modelica_LinearSystems2.DiscreteStateSpace(
          sc,
          dtVar,
          Modelica_LinearSystems2.Types.Method.Trapezoidal);
    (y[:, :, 1],x_continuous[:, :, 1]) :=
      Modelica_LinearSystems2.DiscreteStateSpace.initialResponse(sd, x0, samples);
    else

  for i1 in 1:size(sc.B, 2) loop
        // Loop over inputs

        // time response to plot
    if response == TimeResponse.Impulse then
      u[1, :] := zeros(size(sc.B, 2));
      u[1, i1] := 1;
      sd := Modelica_LinearSystems2.DiscreteStateSpace(
          sc,
          dtVar,
          Modelica_LinearSystems2.Types.Method.ImpulseExact);
    elseif response == TimeResponse.Step then
      u[:, :] := zeros(samples, size(sc.B, 2));
      u[:, i1] := ones(samples);
      sd := Modelica_LinearSystems2.DiscreteStateSpace(
          sc,
          dtVar,
          Modelica_LinearSystems2.Types.Method.StepExact);
    elseif response == TimeResponse.Ramp then
      u[:, :] := zeros(samples, size(sc.B, 2));
      u[:, i1] := 0:dtVar:tSpanVar;
      sd := Modelica_LinearSystems2.DiscreteStateSpace(
          sc,
          dtVar,
          Modelica_LinearSystems2.Types.Method.RampExact);
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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>timeResponse</b>(ss, dt, tSpan, responseType, x0)  </td> </tr>

</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function timeResponse calculates the time responses of a state space system. The type of the time response is defined by the input <b>responseType</b>, i.e. 
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
   Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1],
      B=[1],
      C=[2],
      D=[0]);
  Real Ts=0.1;
  Real tSpan= 0.4;
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  Real x0[1]={0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(ss,Ts,tSpan,response,x0);
//  y[:,1,1]={0, 0.19, 0.3625, 0.518, 0.659}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>


</html> "));
 end timeResponse;

encapsulated function impulseResponse
      "Calculate the impulse time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask2;

algorithm
  (y,t,x_continuous) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
      sc=sc,
      dt=dt,
      tSpan=tSpan,
      response=Modelica_LinearSystems2.Types.TimeResponse.Impulse,
      x0=zeros(size(sc.A, 1)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>impulseResponse</b>(ss, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>impulseResponse</b> calculates the time response of a state space system for impulse imput. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
StateSpace.Analysis.impulseResponse(ss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>



</p>

<h4><font color=\"#008000\">Example</font></h4>
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
  (y,t,x):=StateSpace.Analysis.impulseResponse(ss,Ts,tSpan);
//  y[:,1,1]={2, 1.8097, 1.6375, 1.4816, 1.3406}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 0.9048, 0.8187, 0.7408, 0.6703}
</pre></blockquote>


</html> "));
end impulseResponse;

encapsulated function stepResponse
      "Calculate the step time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask2;

algorithm
  (y,t,x_continuous) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
      sc=sc,
      dt=dt,
      tSpan=tSpan,
      response=Modelica_LinearSystems2.Types.TimeResponse.Step,
      x0=zeros(size(sc.A, 1)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>stepResponse</b>(ss, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>stepResponse</b> calculates the step response of a state space system. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
StateSpace.Analysis.stepResponse(ss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Step, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>


<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
end stepResponse;

encapsulated function rampResponse
      "Calculate the ramp time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask2;

algorithm
  (y,t,x_continuous) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
      sc=sc,
      dt=dt,
      tSpan=tSpan,
      response=Modelica_LinearSystems2.Types.TimeResponse.Ramp,
      x0=zeros(size(sc.A, 1)));

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>rampResponse</b>(ss, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>rampResponse</b> calculates the time response of a state space system for ramp imput u = t. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
StateSpace.Analysis.rampResponse(ss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Ramp, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
end rampResponse;

encapsulated function initialResponse
      "Calculate the time response of a state space system for given initial condition and zero inputs"

      import Modelica;
      import Modelica_LinearSystems2;

  input Real x0[:]=fill(0,0) "Initial state vector";

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask2(redeclare Real y[:,size(sc.C,1),1], redeclare
          Real x_continuous[
                        :,size(sc.A,1),1]);

algorithm
  (y,t,x_continuous) := Modelica_LinearSystems2.StateSpace.Analysis.timeResponse(
      sc=sc,
      dt=dt,
      tSpan=tSpan,
      response=Modelica_LinearSystems2.Types.TimeResponse.Initial,
      x0=x0);

annotation(interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>initialResponse</b>(ss, dt, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initialResponse</b> calculates the time response of a state space system for given initial condition and zero inputs. 
The state space system is transformed to a appropriate discrete state space system and, starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dt.
<blockquote><pre>
StateSpace.Analysis.initialResponse(x0,ss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
StateSpace.Analysis.timeResponse(ss, dt, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.timeResponse\">StateSpace.Analysis.timeResponse</a>
</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
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
      Modelica_LinearSystems2.TransferFunction tf= if StateSpace.Internal.isSISO(ss) then StateSpace.Conversion.toTransferFunction(ss) else TransferFunction(1);

  algorithm
    assert(StateSpace.Internal.isSISO(ss),"System must be SISO but is "+ String(size(ss.B,2)) +"-by-" +String(size(ss.C,1)) +" system");
   result := size(tf.n,1)-1;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>numeratorDegree</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function Analysis.<b>numeratorDegree</b> calculates the degree of the numerator polynomial of the corresponding transfer function. 
The state space system is converted to the transfer function G(s)=N(s)/D(s) with the polynomial N(s) as numerator.
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a> and <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.denominatorDegree\">StateSpace.Analysis.denominatorDegree</a>.
</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
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
    Modelica_LinearSystems2.TransferFunction tf= if StateSpace.Internal.isSISO(ss) then StateSpace.Conversion.toTransferFunction(ss) else TransferFunction(1);

  algorithm
    assert(StateSpace.Internal.isSISO(ss),"System must be SISO but is "+ String(size(ss.B,2)) +"-by-" +String(size(ss.C,1)) +" system");
  result := size(tf.d,1)-1;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>denominatorDegree</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function Analysis.<b>denominatorDegree</b> calculates the degree of the denominator polynomial of the corresponding transfer function. 
The state space system is converted to the transfer function G(s)=N(s)/D(s) with the polynomial D(s) as denominator.
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a> and <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.numeratorDegree\">StateSpace.Analysis.numeratorDegree</a>.
</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end denominatorDegree;

  encapsulated function evaluate
      "Evaluate a the corresponding transfer function at a given (complex) value of s"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;
    input Complex s "Value of s where tf shall be evaluated";
    input Real den_min=0 "|denominator(s)| is limited by den_min";
    output Complex result "= tf(s)";

    protected
    Boolean issiso=StateSpace.Internal.isSISO(ss);
    TransferFunction tf= if issiso then StateSpace.Conversion.toTransferFunction(ss) else TransferFunction(1);
    Complex j=Modelica_LinearSystems2.Math.Complex.j();
    Complex den=Polynomial.evaluateComplex(Polynomial(tf.d), s);
    Real abs_den=Complex.'abs'(den);
  algorithm
    assert(issiso,"System must be SISO but is "+ String(size(ss.B,2)) +"-by-" +String(size(ss.C,1)) +" system");
    den := if abs_den >= den_min then den else -abs_den + 0*j;
    result := Polynomial.evaluateComplex(Polynomial(tf.n), s)/den;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  result </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>evaluate</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the corresponding transfer function of the state space system at a given (complex) value of s.
The state space system is converted to the transfer function G(s)=N(s)/D(s), which is evaluated by calculating the numerator polynomial N(s) and the denominator polynomial D(s).
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a> and <a href=\"Modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end evaluate;

  encapsulated function zerosAndPoles
      "Calculate zeros and poles of the TransferFunction corresponding to a state space representation"
      import Modelica;
      import Modelica_LinearSystems2.Math.Complex;
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

    z:=Polynomial.roots(Polynomial(tf.n));
    p:=Polynomial.roots(Polynomial(tf.d));
    pn:=Polynomial(z);
    pd:=Polynomial(p);
    tf2:=TransferFunction(pn, pd);
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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (z,p,k) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>zerosAndPoles</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function calculates the zeros, poles and gain of the corresponding transfer function of a state space system.
See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">StateSpace.Conversion.toTransferFunction</a> and <a href=\"Modelica://Modelica_LinearSystems2.TransferFunction.Analysis.zerosAndPoles\">TransferFunction.Analysis.zerosAndPoles</a>

</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end zerosAndPoles;

  encapsulated function eigenValues
      "Calculate the eigenvalues of a linear state space system and write them in a complex vector"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss "state space system";
    output Complex eigvalues[size(ss.A, 1)]=Complex.eigenValues(ss.A)
        "eigen values of the system";
  algorithm

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  eigenvalues </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>eigenValues</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Calculate the eigenvalues of a state space system, i.e. the eigenvalues of the system matrix <b>A</b> of a state space system. The output is a complex vector containing the eigenvalues.


</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end eigenValues;

  encapsulated function eigenVectors
      "Calculate the rigth eigenvectors of a linear state space system and write them columnwise in a matrix. Optionally, the eigenvalues are computed"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss "state space system";
    input Boolean onlyEigenvectors=true;
    output Real eigvec[size(ss.A, 1),size(ss.A, 2)]
        "eigen values of the system";
    output Complex eigval[size(ss.A, 1)]=fill(Complex(0), size(ss.A, 1))
        "eigen values of the system";
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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (eigenvectors, eigenvalues) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>eigenVectors</b>(ss, onlyEigenvectors)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Calculate the eigenvectors and optionally (onlyEigenvectors=false) the eigenvalues of a state space system. The output <tt>eigenvectors</tt> is a matrix with the same dimension as matrix <b>ss.A</b>. Just like in <a href=\"Modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>, if the i-th eigenvalue has an imaginary part, then <tt>eigenvectors</tt>[:,i] is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the eigenvector of the i-th eigenvalue.<br>
The eigenvalues are returned as a complex vector <tt>eigenvalues</tt>.


</p>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end eigenVectors;

  encapsulated function invariantZeros
      "Compute invariant zeros of linear state space system"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

    input StateSpace ss "Linear system in state space form";

    output Complex Zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
    Integer n=10;
    Integer m;
    Integer p;
    Real Ar[:,:];
    Real Br[:,:];
    Real Cr[:,:];
    Real Dr[:,:];

    Real Af[:,:];
    Real Bf[:,:];
    Real AfBf[:,:];

    Real V2[:,:];
    Real Vf[:,:];
    Real R[:,:];

    Integer na;
   Real alphaReal[:];
    Real alphaImag[:];
    Real beta[:];
    Integer info;
    Real beta_small=100*Modelica.Constants.eps;
    Real normB=max(Modelica.Math.Matrices.norm(ss.B, p=1),beta_small);
    Real normA=max(Modelica.Math.Matrices.norm(ss.A, p=1),beta_small);

  algorithm
    if min(size(ss.B)) == 0 or min(size(ss.C)) == 0 then
      Zeros := fill(Complex(0), 0);
    else
      (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(ss.A, ss.B, ss.C, ss.D);
      if n > 0 then
        (Ar,Br,Cr,Dr,n,m,p) := StateSpace.Internal.reduceRosenbrock(transpose(Ar), transpose(Cr), transpose(Br), transpose(Dr));
      end if;
      if n == 0 then
        Zeros := fill(Complex(0), 0);
      else
        (,R,,V2) := Matrices.QR(Matrices.fliplr(transpose([Cr,Dr])));
        Vf := Matrices.fliplr(V2);
        AfBf := [Ar,Br]*Vf;
        Af := AfBf[:, 1:size(Ar, 2)];
        Bf := Vf[1:size(Ar, 1), 1:size(Ar, 2)];

        (alphaReal,alphaImag,beta,,,info) := LAPACK.dggev(Af, Bf, n);
        assert(info == 0, "Failed to compute invariant zeros with function invariantZeros(..)");

        Zeros := fill(Complex(0), size(beta, 1));
        normB:=max(Modelica.Math.Matrices.norm(Bf), beta_small);
        normA:=max(Modelica.Math.Matrices.norm(Af, p=1), beta_small);

  // If beta[i] is zero, then zero i is infinite.
        for i in 1:size(beta, 1) loop
          if beta[i] >= normB*1e-3 then
       // finite eigenvalue
            Zeros[i].re := if abs(alphaReal[i]) >= normB*1e-12 then alphaReal[i]/
              beta[i] else 0;
            Zeros[i].im := if abs(alphaImag[i]) >= normB*1e-12 then alphaImag[i]/
              beta[i] else 0;
          end if;
        end for;

      end if;
    end if;
    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  zeros </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>invariantZeros</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
s  that make the Rosenbrock matrix of the system
</p>
<pre>
    | s<b>I-A</b>   <b>-B</b> |
    |           |
    | <b>C</b>       <b>D</b> |

</pre>
singular.
<p>
This function applies the algorithm described in [1] where the system (<b>A</b>, <b>B</b>, <b>C</b>, <b>D</b>) is reduced to a new system (<b>A</b>r, <b>B</b>r <b>C</b>r, <b>D</b>r) with the same zeros and with <b>D</b>r of full rank.
</p>


</p>

<h4><font color=\"#008000\">Example</font></h4>
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

<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center>  Emami-Naeini, A. and Van Dooren, P. </td>  <td> \"Computation of Zeros of Linear Multivariable Systems\"  </td> <td> Automatica, 18, pp. 415-430, 1982. </td></tr>
</table>
</html> "));
  end invariantZeros;

  encapsulated function dcGain
      "Return steady state gain matrix K (for a stable system: K[i,j] = value of y[i] at infinite time for a step input of u[j])"

    import Modelica;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.Math.Complex;
    import Modelica_LinearSystems2.Math.Matrices;
    import Modelica_LinearSystems2.Math.Matrices.LAPACK;

    input StateSpace ss "Linear system in state space form";

    output Real K[size(ss.C,1), size(ss.B,2)]
        "DC gain matrix K (= G(s=0) = D - C*inv(A)*B)";
    output Boolean finite
        "= true, if K is finite; = false, if K is infinite (K=fill(Modelica.Constants.inf,..) returned)";

    protected
    Integer nx = size(ss.A,1);
    Integer nu = size(ss.B,2);
    Integer ny = size(ss.C,1);
    Real X[nx,nu];
    Integer rank;
  algorithm
    finite :=true;
    if nu == 0 or ny == 0 then
      K := fill(0.0, ny, nu);
    else
      (X, rank) := Modelica_LinearSystems2.Math.Matrices.leastSquares2(ss.A, ss.B);
      // Determine whether A*X-B=0 is not fulfilled (since no unique solution)
      if rank < nx then
         if Modelica.Math.Matrices.norm(ss.A*X-ss.B, p=Modelica.Constants.inf)
            >= 1000*Modelica.Constants.eps then
            finite :=false;
         end if;
      end if;

      if finite then
         // A*X - B = 0:
         K :=ss.D - ss.C*X;
      else
         // The least squares solution does not fulfill A*X - B = 0
         K :=fill(Modelica.Constants.inf, ny, nu);
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

<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>

<p>
is solved for <b>y</b> under steady state conditions, i.e.,
</p>

<pre>
   der(<b>x</b>) = <b>0</b>
</pre>

<p>
resulting in
</p>

<pre>
    <b>y</b> = ( <b>D</b> + <b>C</b>*inv(<b>A</b>)*<b>B</b> )*<b>u</b>
      = <b>K</b>*<b>u</b>
</pre>

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



</html> "));
  end dcGain;

  encapsulated function isControllable
      "Check controllability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.StateSpace.Internal;

    input StateSpace ss;
    input Modelica_LinearSystems2.Types.StaircaseMethod method=Modelica_LinearSystems2.Types.StaircaseMethod.SVD;

    output Boolean controllable;
  algorithm

    controllable := if StateSpace.Internal.isSISO(ss) then 
      StateSpace.Internal.isControllableSISO(ss) else StateSpace.Internal.isControllableMIMO(ss,method);

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  controllable </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>isControllable</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function StateSpace.Analysis.<b>isControllable</b> checks the observability of a state space system. Therefore, the system is transformed into staircase form, i.e. the system matrix <b>H</b> of the transformed system has block upper Hessenberg form:
<blockquote><pre>
     | H11    H12     H13    ...     H1k |           
     | H21    H22     H23    ...     H2k |          
 <b>H</b> = |  0     H32     ...    ...     ... |
     | ...    ...     ...    ...     ... |            
     |  0     ...      0    Hk,k-1   Hkk |             
</pre>
</blockquote>
where, if <b>H</b>k,k-1 has full rank, indicating whether the system is controllable or not.<br>


For single input systems the staircase form is a usual upper Hessenberg form, i.e. th blocks are of dimension one.<br>
The boolean input <b>method</b> defines for multi output systems the method to generate the staircase form of the system, whereas Types.StaircaseMethod.QR and Types.StaircaseMethod.SVD denotes QR-factorization and singular value decomposition respectively. Since staircase algorithm contains rank decisions QR-factorization should be restricted to well conditioned systems of lower order (<5). Default is SVD.<br>

Since controllability is dual to observability of the dual system (A', C', B', D'), proof of <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.isObservable\">observability</a> is referred to proof of controllability of the dual system.<br>



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
  end isControllable;

  encapsulated function isObservable
      "Check observability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;
    input Modelica_LinearSystems2.Types.StaircaseMethod method=
        Modelica_LinearSystems2.Types.StaircaseMethod.SVD;

    output Boolean observable;
  algorithm

    observable := if StateSpace.Internal.isSISO(ss) then 
      StateSpace.Internal.isObservableSISO(ss) else 
      StateSpace.Internal.isObservableMIMO(ss, method);

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  observable </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>isObservable</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function StateSpace.Analysis.<b>isObservable</b> checks the observability of a state space system. Since observability is dual to controllability of the dual system (A', C', B', D'), proof of observability is referred to proof of <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.isControllable\">controllability</a> of the dual system.<br>
The boolean input <b>method</b> defines for multi output systems the method to generate the staircase form of the system, whereas Types.StaircaseMethod.QR and Types.StaircaseMethod.SVD denotes QR-factorization and singular value decomposition respectively. Since staircase algorithm contains rank decisions QR-factorization should be restricted to  well conditioned systems of lower order (<5). Default is SVD.<br>


</p>

<h4><font color=\"#008000\">Example</font></h4>
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
</html> "));
  end isObservable;

  encapsulated function isStabilizable
      "Check stabilizability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

    output Boolean stabilizable;

  algorithm
    if StateSpace.Internal.isSISO(ss) then
      stabilizable := StateSpace.Internal.isStabilizableSISO(ss);
    else
      stabilizable := StateSpace.Internal.isStabilizableMIMO(ss);
    end if;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  stabilizable </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>isStabilizable</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function checks whether a state space system is stabilizable or not.<br>
A system is stabilizable for the continuous-time case if all of the uncontrollable eigenvalues have negative real part.
Therefore, a controllable system is always stabilizable.
<br>
To check stabilizability, staircase algorithm is used to separate the controllable subspace from the uncontrollable subspace.
Then, the uncontrollable poles are checked to be stable, i.e. to have negative real parts.

</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
  end isStabilizable;

  encapsulated function isDetectable
      "Check detectability of a state space system"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

    output Boolean detectable;

  algorithm
   if StateSpace.Internal.isSISO(ss) then
    detectable := StateSpace.Internal.isDetectableSISO(ss);
    else
     detectable := StateSpace.Internal.isDetectableMIMO(ss);
     end if;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  detectable </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>isDetectable</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function checks whether a state space system is detectable or not.<br>
A system is detectable for the continuous-time case if all of the unobservable eigenvalues have negative real part.
Therefore, a observable system is always detectable.
<br>
To check detectability, staircase algorithm is used to separate the observable subspace from the unobservable subspace.
Then, the unobservable poles are checked to be stable, i.e. to have negative real parts.




</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
  end isDetectable;

  encapsulated function controllabilityMatrix
      "Compute the controllability matrix [B, A*B, ..., A^(n-1)*B] of a state space system"

      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;
    output Real cm[size(ss.B, 1),size(ss.A, 2)*size(ss.B, 2)];

  algorithm
    if size(ss.A, 2) == 0 then
      cm := fill(
          0,
          size(ss.B, 1),
          0);
    else
      cm[:, 1:size(ss.B, 2)] := ss.B;

      for i in 2:size(ss.A, 1) loop
        cm[:, ((i - 1)*size(ss.B, 2) + 1):(i*size(ss.B, 2))] := ss.A*cm[:, ((i
           - 2)*size(ss.B, 2) + 1):((i - 1)*size(ss.B, 2))];
      end for;
    end if;

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  Q </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>controllabilityMatrix</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function calculates the controllability matrix
<blockquote><pre>
  <b>Q</b> = [<b>B</b>, <b>A</b>*<b>B</b>, ..., <b>A</b>^(n-1)*<b>B</b>]
</pre>
</blockquote>

of the system
<blockquote><pre>
  der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
      <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;

</pre>
</blockquote>



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
  end controllabilityMatrix;

  encapsulated function observabilityMatrix
      "Compute the observability matrix of a state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;
    output Real om[size(ss.A, 1)*size(ss.C, 1),size(ss.C, 2)];

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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  Q </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>observabilityMatrix</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function calculates the observability matrix
<blockquote><pre>
   <b>Q</b> = [<b>C</b>; <b>C</b>*<b>A</b>; ...; <b>C</b>*<b>A</b>^(n-1)] 
</pre>
</blockquote>

of the system
<blockquote><pre>
  der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>;
      <b>y</b>  = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>;

</pre>
</blockquote>



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
  end observabilityMatrix;

end Analysis;

  encapsulated package Design "Functions for state space controller design"
    encapsulated function assignPolesSI
      "Pole placement for single input systems using Ackermann's formula."

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      input Modelica_LinearSystems2.Math.Complex p[size(ss.A, 1)]
        "Vector of desired poles";
      output Real k[size(ss.A, 1)] "Feedback gain matrix";

    protected
      Real cm[size(ss.B, 1),size(ss.A, 1)*size(ss.B, 2)];
      Modelica_LinearSystems2.Math.Polynomial poly;
      Real Y[size(ss.A, 1),size(ss.A, 2)];
      Real X[:,:];
      Modelica_LinearSystems2.Math.Complex p_actual[size(p, 1)];
      Modelica_LinearSystems2.Math.Complex p_sorted[size(p, 1)];
      Real poleError;
      Modelica_LinearSystems2.Math.Complex smaller;
    algorithm
      assert(size(ss.B,2)==1,"System must be SI but has "+ String(size(ss.B,2)) +" inputs");
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
          if Modelica_LinearSystems2.Math.Complex.'abs'(p_sorted[i1]) >
              Modelica_LinearSystems2.Math.Complex.'abs'(p_sorted[i2]) then
            smaller := p_sorted[i2];
            p_sorted[i2] := p_sorted[i1];
            p_sorted[i1] := smaller;
          end if;
        end for;
      end for;

      p_actual := Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A - ss.B*
        transpose(matrix(k)));
       // sort actual eigenvalues
      for i1 in 1:size(p_actual, 1) loop
        for i2 in (1 + i1):size(p_actual, 1) loop
          if Modelica_LinearSystems2.Math.Complex.'abs'(p_actual[i1]) >
              Modelica_LinearSystems2.Math.Complex.'abs'(p_actual[i2]) then
            smaller := p_actual[i2];
            p_actual[i2] := p_actual[i1];
            p_actual[i1] := smaller;
          end if;
        end for;
      end for;

      // check for poles that have an error of more than 10%
      for i in 1:size(p_sorted, 1) loop
        if (Modelica_LinearSystems2.Math.Complex.'abs'(p_sorted[i]) <> 0) then
           poleError := Modelica_LinearSystems2.Math.Complex.'abs'((p_sorted[i] -
            p_actual[i]))/Modelica_LinearSystems2.Math.Complex.'abs'(p_sorted[i]);

          if poleError > 0.1 then
              Modelica.Utilities.Streams.print("Warning: Pole location of pole " +
              String(p_sorted[i]) + " has an error of " + String(100*poleError)
               + "%. (Is " + String(p_actual[i]) + ")");

          end if;
        end if;
      end for;

    end assignPolesSI;

    encapsulated function assignPolesMI
      "Pole assigment design algorithm for multi input systems"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Math.Matrices;

      input StateSpace ss "state space system";

      input Complex gamma[:]=fill(Complex(0), 0) "Designed Poles";
    //  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
      input Real alpha=-1e10
        "maximum admissible value for real parts(continuous) or for moduli (discrete) of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
      input Real tolerance=Modelica.Math.Matrices.norm(ss.A, 1)*1e-12
        "The tolerance to be used in determining the controllability of (A,B)";
      input Boolean calculateEigenvectors=false
        "Calculate the eigenvectors X of the closed loop system when true";

      output Real K[size(ss.B, 2),size(ss.A, 1)]
        "State feedback matrix assigning the desired poles";
      output Real S[:,:] "Closed loop System matrix";
      output Complex po[size(ss.A, 1)] "poles of the closed loop system";
      output Integer nfp
        "number of eigenvalues that are not modified with respect to alpha";
      output Integer nap "number of assigned eigenvalues";
      output Integer nup "number of uncontrollable eigenvalues";
      output Complex X[size(ss.A, 1),size(ss.A, 1)]
        "eigenvectors of the closed loop system";

    protected
      Real A_rsf[size(ss.A, 1),size(ss.A, 2)];
      Real B_rsf[size(ss.B, 1),size(ss.B, 2)];
      Real Q[size(ss.A, 1),size(ss.A, 1)];
      Real Ks1[:,:];
      Real Ks2[:,:];
      Real Q2[:,:];
      Real A_rsf_1[:,:];
      Real Q1[:,:];
      Boolean select[:];
      Boolean rselectA[:];
      Real Z[:,:] "orthogonal transformation matrix";
      Real ZT[:,:] "orthogonal transformation matrix";
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
      Integer ng=size(gamma,1);
      Integer nr "Differenz between rpA and rpg; Sign(rpA-rpg)*(rpA-rpg)";

      Real alphaReal[size(ss.A, 1)]
        "Real part of eigenvalue=alphaReal+i*alphaImag";
      Real alphaImag[size(ss.A, 1)]
        "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

      Complex SS[:,:];
      Complex Xj[:,:];
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
      (A_rsf,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSF2(
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
        rpA := n-nccA;
        nccA := div(nccA, 2);

      // reorder gamma and A_rsf
      (gammaReordered,rpg) := Modelica_LinearSystems2.Internal.reorderZeros(gamma);
      gammaReordered := Complex.Vectors.reverse(gammaReordered);
      nccg := div(size(gammaReordered, 1) - rpg, 2);
      ncc := min(nccA, nccg);
      rp := min(rpA, rpg);
      if nccA > 0 then
        (A_rsf[nfp + 1:n, nfp + 1:n],Q2) := Matrices.LAPACK.dtrsen(
            "E",
            "V",
            rselectA,
            A_rsf[nfp + 1:n, nfp + 1:n],
            identity(n - nfp));//The Schur vector matrix is identity, since A_rsf already has Schur form

        A_rsf[1:nfp, nfp + 1:n] := A_rsf[1:nfp, nfp + 1:n]*Q2;
        B_rsf[nfp + 1:n, :] := transpose(Q2)*B_rsf[nfp + 1:n, :];
        ZT[nfp + 1:n, :] := transpose(Q2)*ZT[nfp + 1:n, :];
      end if;

      // main algorithm
      K := zeros(size(ss.B, 2), size(ss.A, 1));
      counter := nfp + 1;
      counter2 := 1;

      for i in 1:rp loop // 1x1 blocks; real system pole and real assigned poles; take the next eigenvalue in the
                           // diagonal of the Schur form and search the nearest pole in the set of the real poles to assign
          dist := Modelica.Constants.inf;
          for ii in i:rpg loop // looking for nearest pole and reorder gamma
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
          K := K + [zeros(size(Ks1, 1), size(K, 2) - 1),Ks1]*ZT;
          A_rsf := A_rsf - B_rsf*[zeros(size(Ks1, 1), size(K, 2) - 1),Ks1];
          select := fill(false, n - counter + 1);
          select[n - counter + 1] := true;

          (A_rsf[counter:n, counter:n],Q1) := Matrices.LAPACK.dtrsen(
            "E",
            "V",
            select,
            A_rsf[counter:n, counter:n],
            identity(n - counter + 1));//The Schur vector matrix is identity, since A_rsf already has Schur form

          A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q1;
          B_rsf[counter:n, :] := transpose(Q1)*B_rsf[counter:n, :];
          ZT[counter:n, :] := transpose(Q1)*ZT[counter:n, :];
          counter := counter + 1;
          counter2 := counter2 + 1;
        end for;

        if counter2<rpg and counter2>rpA then //System has less real eigenvalues than real assigned poles
        for i in 1:div(rpg - rpA, 2) loop // 2x2 blocks; complex pair of system poles and 2 real assigned poles; take the next complex pair
                                          // (Schur bump) in the diagonal of the Schur form and search the two nearest poles in the set of the
                                          // remaining real assigned poles
          dist := Modelica.Constants.inf;
          evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);//positive imaginary part of the complex system pole pair
          for ii in 2*(i - 1) + 1:2:rpg - rpA loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag;
            end if;
          end for;
          h := gammaReordered[ng - rp - 2*(i - 1)];
          gammaReordered[ng - rp - 2*(i - 1)] := gammaReordered[iii];
          gammaReordered[iii] := h;
          dist := Modelica.Constants.inf;
          for ii in 2*(i - 1) + 1:2:rpg - rpA loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + evImag;
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

          K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
          A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
          select := fill(false, n - counter + 1);
          select[n - counter:n - counter + 1] := {true,true};

          (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
            "E",
            "V",
            select,
            A_rsf[counter:n, counter:n],
            identity(n - counter + 1)); //The Schur vector matrix is identity, since A_rsf already has Schur form

          A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
          B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
          ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
          counter := counter + 2;
          counter2 := counter2 + 2;
        end for;
      end if;

      if counter2>rpg and counter2<rpA then//System has more real eigenvalues than real assigned poles
        for i in 1:div(rpA - rpg, 2) loop// 2x2 blocks; 2 real system poles and a pair of complex assigned poles; take the next two real
                                          // eigenvalues in the diagonal of the Schur form and search the complex pole pair of the assigned poles
                                          // which is nearest to the two real poles
          dist := Modelica.Constants.inf;
          for ii in 2*(i - 1)+1:2:rpA - rpg loop
            if abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[n - 1, n - 1] -
              gammaReordered[ng - rp - ii + 1].re) + abs(gammaReordered[ng - rp - ii + 1].im) < dist then
              iii := ng - rp - ii + 1;
              dist := abs(A_rsf[n, n] - gammaReordered[ng - rp - ii + 1].re)
                 + abs(gammaReordered[ng - rp - ii + 1].im) + abs(A_rsf[
                n - 1, n - 1] - gammaReordered[ng - rp - ii + 1].re) +
                abs(gammaReordered[ng - rp - ii + 1].im);
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

          K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
          A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
          select := fill(false, n - counter + 1);
          select[n - counter:n - counter + 1] := {true,true};

          (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
            "E",
            "V",
            select,
            A_rsf[counter:n, counter:n],
            identity(n - counter + 1)); //The Schur vector matrix is identity, since A_rsf already has Schur form

          A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
          B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
          ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
          counter := counter + 2;
          counter2 := counter2 + 2;
          Modelica.Utilities.Streams.print("counter2Case3 = " + String(counter2));
        end for;
      end if;

      for i in 1:ncc loop // 2x2 blocks; 2 complex system poles and two complex assigned poles; take the next complex
                          // system pole pair (next Schur bump) in the diagonal of the Schur form and search the complex
                          //  assigned pole pair which is nearest
        dist := Modelica.Constants.inf;
        evImag := sqrt(-A_rsf[n - 1, n]*A_rsf[n, n - 1]);//positive imaginary part of the complex system pole pair
        for ii in 2*(i - 1) + 1:2:2*ncc loop
          if abs(A_rsf[n, n] - gammaReordered[2*ncc - ii + 1].re) + abs(evImag -
              abs(gammaReordered[2*ncc - ii + 1].im)) < dist then
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
        K := K + [zeros(size(Ks2, 1), size(K, 2) - 2),Ks2]*ZT;
        A_rsf := A_rsf - B_rsf*[zeros(size(Ks2, 1), size(K, 2) - 2),Ks2];
        select := fill(false, n - counter + 1);
        select[n - counter:n - counter + 1] := {true,true};

        (A_rsf[counter:n, counter:n],Q2) := Matrices.LAPACK.dtrsen(
          "E",
          "V",
          select,
          A_rsf[counter:n, counter:n],
          identity(n - counter + 1));   //The Schur vector matrix is identity, since A_rsf already has Schur form

        A_rsf[1:counter - 1, counter:n] := A_rsf[1:counter - 1, counter:n]*Q2;
        B_rsf[counter:n, :] := transpose(Q2)*B_rsf[counter:n, :];
        ZT[counter:n, :] := transpose(Q2)*ZT[counter:n, :];
        counter := counter + 2;
        counter2 := counter2 + 2;
      end for;

      S := ss.A - ss.B*K;
      po := Complex.eigenValues(S);

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
        X := Complex.eigenVectors(S);
    //      Modelica_LinearSystems2.Math.Complex.Matrices.print(X,6,"X2");

      end if;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (K, S, po, nfp, nap, nup) </td><td align=center> =  </td>  <td> StateSpace.Design.<b>assignPolesMI</b>(ss, gamma, np, tol, calculateEigenvectors)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
The purpose of this function is to determine the state feedback matrix <b>K</b> for a
given time invariant multi input state system (<b>A</b>,<b>B</b>) such that the
closed-loop state matrix <b>A</b>-<b>B</b>*<b>K</b> has specified eigenvalues. The
feedback matrix <b>K</b> is calculated by factorization following [1]. The algorithm
modifies the eigenvalues sequentially and also allows partial eigenvalue assignment.<br>
<br>


 At the beginning of the algorithm, the feedback matrix <b>K</b> is set to zero (<b>K</b> = <b>0</b>) and the matrix <b>A</b> is
 reduced to an ordered real Schur form by separating its spectrum in two parts

<blockquote><pre>
              | <b>F</b>1  <b>F</b>3|
 <b>F</b> = <b>Q</b>*<b>A</b>*<b>Q</b>' = |       |
              | <b>0</b>   <b>F</b>2|
 </pre>
</blockquote> in such a way, that <b>F</b>1 contains the eigenvalues that will be
retained and <b>F</b>3 contains the eigenvalues going to be modified. On the suggestion
of [1] the eigenvalues <i>evr</i> to be retained are chosen as
 <blockquote><pre>
  evr = {s in C: Re(s) &lt -alpha, alpha &gt =0}
 </pre> </blockquote>
but other specification are conceivable of course.<br>
<br>

Let
 <blockquote><pre>
  <b>G</b> = [<b>G</b>1;<b>G</b>2] = <b>Q</b>*<b>B</b>
 </pre> </blockquote>
with an appropriate partition according to <b>F</b>2. (<b>F</b>2, <b>G</b>2) has to be
controllable.<br>

If the feedback matrix <b>K</b> is taken in a form <blockquote><pre> <b>K</b> = [0, <b>K</b>2]
</pre></blockquote> the special structure of <b>F</b> and <b>K</b> results in a closed loop state
matrix <blockquote><pre>
          |<b>F</b>1 <b>F</b>3 - <b>G</b>1*<b>K</b>2|
<b>F</b> - <b>G</b>*<b>K</b> = |             |
          |0  <b>F</b>2 - <b>G</b>2*<b>K</b>2|
</pre></blockquote> with only the eigenvalues of <b>F</b>2 are modified. This approach to modify
separated eigenvalues is used to sequentially shift one real eigenvalue ore two
complex conjugated eigenvalues stepwise until all assigned eigenvalues are placed.
Therefore, at each step i always the (two) lower right eigenvalue(s) are modified by an
appropriate feedback matrix <b>K</b>i. The matrix <b>F</b> - <b>G</b>*<b>K</b>i remains in real Schur form. The
assigned eigenvalue(s) is (are) then moved to another diagonal position of the real Schur
form using reordering techniques <b>F</b> &lt -- <b>Q</b>i*<b>F</b>*<b>Q</b>i'  and a new block is transferred to the
lower right diagonal position. The transformations are accumulated in <b>Q</b>i and are also
applicated to the matrices <blockquote><pre> <b>G</b> &lt - <b>Q</b>i*<b>G</b> <b>Q</b> &lt - <b>Q</b>i*<b>Q</b> </pre></blockquote>
The eigenvalue(s) to be assigned at  each step is (are) chosen such that the norm of each <b>K</b>i is minimized [1].
<p>



</p>

<h4><font color=\"#008000\">Example</font></h4>
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


<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center>  Varga A.  </td>  <td> \"A Schur method for pole assignment\"  </td> <td> IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519, 1981 </td></tr>
</table>

</html> "));
    end assignPolesMI;

  function kalmanFilter "Design of a Kalman estimator matrix"
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss "Time-continuous system in state space form";
    input Real Q[size(ss.A, 1),size(ss.A, 1)]
        "Covariance Matrix of state noise (n x n), n number of states";
    input Real R[size(ss.C, 1),size(ss.C, 1)]
        "Covariance Matrix of output noise (m x m), m number of inputs";

    public
    output Real L[:,:] "Kalman filter matrix";
    output StateSpace kss(
      redeclare Real A[size(ss.A, 1),size(ss.A, 1)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2) + size(ss.C, 1)],
      redeclare Real C[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real D[size(ss.A, 1),size(ss.B, 2) + size(ss.C, 1)])
        "kalman system";

    protected
    Real AR[:,:]=transpose(ss.A);
    Real BR[:,:]=transpose(ss.C);
    Real CR[:,:]=zeros(1, size(AR, 1));
    Real DR[:,:]=zeros(1, size(BR, 2));
    StateSpace rss=StateSpace(AR, BR, CR, DR)
        "System to calculate the Kalman estimator with lqr algorithm";

    Real AK[size(ss.A, 1),size(ss.A, 1)];
    Real BK[size(ss.B, 1),size(ss.B, 2) + size(ss.C, 1)];
    Real CK[size(ss.A, 1),size(ss.A, 2)];
    Real DK[size(ss.A, 1),size(ss.B, 2) + size(ss.C, 1)];// matrices of the kalman system kss

  algorithm
    (L) := StateSpace.Design.lqr(rss, Q, R, true);
    L := transpose(L);

    AK := ss.A - L*ss.C;
    BK[:, 1:size(ss.B, 2)] := ss.B - L*ss.D;
    BK[:, size(ss.B, 2) + 1:size(BK, 2)] := L;
    CK := identity(size(ss.A, 1));
    DK := zeros(size(ss.A, 1),size(ss.B, 2) + size(ss.C, 1));

    kss := StateSpace(AK, BK, CK, DK);

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (L, kss) </td><td align=center> =  </td>  <td> StateSpace.Design.<b>kalmanFilter</b>(ss, Q, R)  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This functions designs the kalman-bucy filter, that reconstructs
plant states and system output without noise.
As input it uses the plant input and output.
</p>
<p>
Noise affects the plant states via q(t)
</p>
<blockquote>dx/dt = Ax + Bu + q(t)
</blockquote>
<p>
The plant output is affected by r(t)
<blockquote>y = Cx + Du + r(t)
</blockquote>
<p>
The covariance matrices of q and r have to be given via Q and R, respectively.
 
<p>
The filter uses an observer that tries to reconstruct the original behaviour. Its states and outputs are trailed with a hat (^)<br>
The observer is controlled by feedback of the output difference y - y^ (y^= Cx^+ Du)
over a Matrix L, such that x^
</p>
<blockquote>dx^/dt = (A - LC) x^ + (B - LD)u + Ly
</blockquote>
<p> 
follows the plant state x. L is designed to minimize noise in states and inputs.
L is calculated from a Riccati Equation
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
Since the controller approach was made to provide the estimated states, the representation of the ooutput kss is such that
<blockquote><pre>
   y^ = x^ 
</pre></blockquote>
i.e., kss:
<blockquote><pre>
  .                            |u|
  x^  = [A-LC] x^ + [B-LD , L] | |
                               |y|
              |u|
  y^ = Ix^ + 0| |
              |y|
            
</pre></blockquote>
i.e.
<blockquote><pre>
  C^ = I,   D^ = 0   with appropriate sizes, i.e. size(C^) = {nx,nx},  size(D^) = {nx, nu+ny}
</pre></blockquote>

<p>
Since the calculation of a Kalman filter is the dual problem of lqr calculation function
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Design.lqr\">Modelica_LinearSystems2.StateSpace.Design.lqr</a>
is used to solve the Riccati euation.<br>
The algebraic Riccati equation is solved by using the Schur algorithm 
<a href=\"Modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a>.
</p>
 
<h4><font color=\"#008000\">Example</font></h4>
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
 
 
 
 
 
</html>"),   DymolaStoredErrors);
  end kalmanFilter;

    encapsulated function lqr "LQR design algorithm"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math;

      input StateSpace ss "open loop system";
      input Real Q[size(ss.A, 1),size(ss.A, 2)]=identity(size(ss.A, 1))
        " state weighting matrix";
      input Real R[size(ss.B, 2),size(ss.B, 2)]=identity(size(ss.B, 2))
        " input weighting matrix";
    protected
      input Boolean iscontinuousSystem=true;
    public
      output Real K[size(ss.B, 2),size(ss.A, 1)] "Feedback gain matrix";
      output StateSpace sslqr(
        redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1),size(ss.D, 2)]) "closed loop system";

      output Real S[size(ss.A, 1),size(ss.A, 1)]
        "solution of the Riccati equation";

      output Math.Complex ev[:];

    algorithm
      if min(size(ss.A,1),size(ss.A,2))>0 then
      assert(StateSpace.Analysis.isControllable(ss),"System in function \"Modelica_LinearSystems2.StateSpace.Design.lqr\" has to be controllable");
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
        K := Math.Matrices.solve2(R + transpose(ss.B)*S*ss.B, transpose(ss.B)*S*
          ss.A);
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
        ev := fill(Modelica_LinearSystems2.Math.Complex(0), 0);
      end if;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (K, sslqr, X, ev) </td><td align=center> =  </td>  <td> StateSpace.<b>lqr</b>(ss, Q, R, true)  </td> </tr>

</table>
<h4><font color=\"#008000\">Description</font></h4>
The optimal and stabilizing gain matrix <b>K</b> for a state-feedback law <b>u</b> = -<b>K</b>*<b>x</b>
is designed such that the cost function
<p>
 <blockquote><pre>
J = Integral {<b>x</b>'*<b>Q</b>*<b>x</b> + <b>u</b>'*<b>R</b>*<b>u</b>} dt
</pre></blockquote>
of the continuous time case or
<p>
 <blockquote><pre>
Jd = Sum {<b>x</b>'k*<b>Q</b>*<b>x</b>k + <b>u</b>'k*<b>R</b>*<b>u</b>k}
</pre></blockquote>
<p>
of the discrete time case is minimized. The cases are chosen by the input <b>iscontinuousSystem</b> This is done by solving
the continuous-time algebraic Riccati equation (CARE)
<p>
<blockquote><pre>
                       -1
 <b>Q</b> + <b>A</b>'*<b>X</b> + <b>X</b>*<b>A</b> - <b>X</b>*<b>B</b>*<b>R </b>*<b>B</b>'*<b>X</b> = <b>0</b>
</pre></blockquote>
or the discrete-time algebraic Riccati equation (DARE)
<blockquote><pre>
                                 -1
 <b>X</b> - <b>A</b>'*<b>X</b>*<b>A</b> + <b>A</b>'*<b>X</b>*<b>B</b>*(<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)  *<b>B</b>'*<b>X</b>*<b>A</b> - <b>Q</b> = <b>0</b>
</pre>
</blockquote>
<p>
for <b>X</b> using the Schur vector approach. See <a href=\"Modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a> and <a href=\"Modelica://Modelica_LinearSystems2.Math.Matrices.dare\">dare</a> respectively for more details.
<p>
The gain matrix <b>K</b> of the continuous-time case is calculated from
<p>
<blockquote><pre>
     -1
<b>K</b> = <b>R</b>  *<b>B</b>'*<b>X</b>
</pre></blockquote>
or from
<blockquote><pre>
               -1
<b>K</b> = (<b>R</b> + <b>B</b>'*<b>X</b>*<b>B</b>)  <b>B</b>'*<b>X</b>*<b>A</b>
</pre></blockquote>
for the discrete-time case.
<p>
The output state space system sslqr represents the closed loop system
<blockquote><pre>
  .
  <b>x</b> = [<b>A</b> - <b>BK</b>] <b>x</b> + <b>Bu</b>

  <b>y</b> = [<b>C</b> - <b>DK</b>] <b>x</b> + <b>Du</b>

</pre></blockquote>
<p>
The output S is the solution of the Riccati equation

<p>
The eigenvalues of the closed loop system <b>A</b> - <b>B</b>*<b>K</b> are computed as complex output ev.
</p>
<h4><font color=\"#008000\">Example</font></h4>
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

</html> "));
    end lqr;

    encapsulated function lqg "LQG design algorithm"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math;

      input StateSpace ss "open loop system";
      input Real Q[size(ss.A, 1),size(ss.A, 2)]=identity(size(ss.A, 1))
        " state weighting matrix";
      input Real R[size(ss.B, 2),size(ss.B, 2)]=identity(size(ss.B, 2))
        " input weighting matrix";
      input Real V[size(ss.C, 1),size(ss.C, 1)]=identity(size(ss.C, 1))
        " covariance output noise matrix";
      input Real W[size(ss.A, 1),size(ss.A, 1)]=identity(size(ss.A, 1))
        " covariance input noise matrix";

      input Boolean iscontinuousSystem=true;

      output Real Kc[size(ss.B, 2),size(ss.A, 1)]
        "Controller feedback gain matrix";
      output Real Kf[size(ss.A, 1),size(ss.C, 1)] "Kalman feedback gain matrix";

      output StateSpace sslqg(
        redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1),size(ss.D, 2)]) "closed loop system";

    protected
     Real AR[:,:]=transpose(ss.A);
     Real BR[:,:]=transpose(ss.C);
     Real CR[:,:]=zeros(1, size(AR, 1));
     Real DR[:,:]=zeros(1, size(BR, 2));  // System for lqr
     StateSpace rss=StateSpace(
         AR,
         BR,
         CR,
         DR);
      Real Sc[size(ss.A, 1),size(ss.A, 1)]
        "solution of the Riccati equation, controller";
      Real Sf[size(ss.A, 1),size(ss.A, 1)]
        "solution of the Riccati equation, filter";

    algorithm
      if min(size(ss.A,1),size(ss.A,2))>0 then
      assert(StateSpace.Analysis.isControllable(ss),"System in function \"Modelica_LinearSystems2.StateSpace.Design.lqg\" has to be controllable");

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
        Kc := Math.Matrices.solve2(R + transpose(ss.B)*Sc*ss.B, transpose(ss.B)*Sc*
          ss.A);
      end if;

      assert(StateSpace.Analysis.isObservable(ss),"System in function \"Modelica_LinearSystems2.StateSpace.Design.lqg\" has to be observable");
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
      sslqg.D := zeros(size(ss.D,1),size(ss.D,2));

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (Kc, Kf, sslqg) </td><td align=center> =  </td>  <td> StateSpace.<b>lqg</b>(ss, Q, R, V, W)  </td> </tr>
</table>

<h4><font color=\"#008000\">Description</font></h4>
<p>
This function calculates matrices <b>K</b>c and <b>K</b>f for linear quadratic gaussian problem (LQG), i.e. the minimization of the expected value of a cost function in consideration of stochastically disturbed states and outputs of the system
</p>
<blockquote><pre>
d<b>x</b>/dt = <b>A</b><b>x</b> + <b>B</b><b>u</b> + <b>w</b>
<b>y</b> = <b>C</b><b>x</b> + <b>D</b><b>u</b> + <b>v</b>
</pre></blockquote>
<p>
The noise <b>w</b>(t) and <b>v</b>(t) are supposed to be both white, Gaussian zero-mean, stationary stochastic processes with positive semidefinte covariance matrix <b>W</b>
<p>
<blockquote>
E[<b>w</b>(t)*<b>w</b>'(tau)] = <b>W</b>*delta(t-tau)
</blockquote>
<p>
and positive covariance matrix <b>V</b>
<p>
<blockquote>
E[<b>v</b>(t)*<b>v</b>'(tau)] = <b>V</b>*delta(t-tau).
</blockquote>
<p>
E[s] denotes the expected value of a signal s.
<p>
The LQG approach combines the deterministic <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Design.lqr\">LQR</a> approach and <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Design.kalmanFilter\">Kalman filter</a> principle to estimate stochastically disturbed systems, such that input <b>u</b>(t) is given by
<p>
<blockquote>
<b>u</b>(t) = -<b>K</b>c<b>x</b>^(t)
</blockquote>
<p>
where <b>K</b>c is a lqr feedback matrix and x^(t) the reconstructed state vector estimated by a Kalman filter.
<p>
Since, the considered problem is stochastic, the objective function to minimize is an expected value
<p>
<blockquote><pre>
           1      T
J = lim   ---- E[Integral (<b>x</b>'*<b>Q</b>*<b>x</b> + <b>u</b>'*<b>R</b>*<b>u</b>)dt],
  (T->inf) 2T    -T
</pre></blockquote>
<p>
where the weighting matrices <b>Q</b> and <b>R</b> are, respectively, symmetric positive semidefinite and positive definite.
<p>
The feedback matrix Kc is calculated by
<blockquote><pre>
      -1
<b>K</b>c = <b>R</b> *<b>B</b>'*<b>X</b>c,
</pre></blockquote>
<p>
where <b>X</b>c satisfying the continuous-time algebraic Riccati equation (<a href=\"Modelica://Modelica_LinearSystems2.Math.Matrices.care\">care</a>)
<blockquote><pre>
                          -1
 <b>Q</b> + <b>A</b>'*<b>X</b>c + <b>X</b>c*<b>A</b> - <b>X</b>c*<b>B</b>*<b>R </b>*<b>B</b>'*<b>X</b>c = <b>0</b>.
</pre></blockquote>
<p>
The matrix <b>K</b>f of the filter problem to generate the estimated state vector <b>x</b>^(t) is given by
<p>
<blockquote>
<b>K</b>f = <b>X</b>f*<b>C</b>T*<b>V</b>-1,
</blockquote>
<p>
where <b>X</b>f is satisfying the continuous-time algebraic Riccati equation
<blockquote><pre>
                           -1
 <b>W</b> + <b>A</b>*<b>X</b>f + <b>X</b>f*<b>A</b>' - <b>X</b>f*<b>C</b>'*<b>V </b>*<b>C</b>*<b>X</b>f = <b>0</b>.
</pre></blockquote>

The vector <b>x</b>^(t) satisfies the differential equation
<p>
<blockquote><pre>
.
<b>x</b>^(t) = (<b>A</b> - <b>K</b>f<b>C</b>)<b>x</b>^(t) + (<b>B</b> - <b>K</b>f<b>D</b>)<b>u</b>(t) + <b>K</b>f<b>y</b>(t)
</pre></blockquote>
<p>
Combining the equation state feedback and state estimation, the state vector <b>x</b>(t) and the estimated state vector <b>x</b>^(t) are given by
<blockquote><pre>
  .
 |<b>x</b> |   | <b>A</b>         -<b>B</b><b>K</b>c      |  |<b>x</b> |   | <b>I</b>   <b>0</b> |  | <b>w</b> |
 |  | = |                     |  |  | + |       |  |   |
 |<b>x</b>^|   | <b>K</b>f<b>C</b>   <b>A</b> - <b>B</b><b>K</b>c - <b>K</b>f<b>C</b> |  |<b>x</b>^|   | <b>0</b>  <b>K</b>f |  | <b>v</b> |.

</pre></blockquote>
<br>
Finally, the output sslqg represents the estimated system with <b>y</b>(t), the output of the real system, as the input
<blockquote><pre>
 .
 <b>x</b>^ = [<b>A</b> - <b>K</b>f<b>C</b> - <b>B</b><b>K</b>c + <b>K</b>f<b>D</b><b>K</b>c]*<b>x</b>^ + <b>K</b>f*<b>y</b>

 <b>y</b>^ = [<b>C</b> - <b>D</b><b>K</b>c] <b>x</b>^

</pre></blockquote>

</p>
<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  StateSpace ss=StateSpace(
      A=[-0.02, 0.005, 2.4,  -32; -0.14,  0.44,  -1.3,  -30; 0,  0.018,  -1.6,  1.2; 0, 0, 1, 0],
      B=[0.14,  -0.12; 0.36, -8.6; 0.35, 0.009; 0, 0],
      C=[0, 1, 0, 0; 0, 0, 0, 57.3],
      D=[0,0; 0,0]);
 
   Real Q[:,:] = transpose(ss.C)*ss.C \" state weighting matrix\";
   Real R[:,:] = identity(2) \" input weighting matrix\";
   Real V[:,:] = identity(2) \" covariance output noise matrix\";
   Real W[:,:] = ss.B*transpose(ss.B) \" covariance input noise matrix\";
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

</html> "));
    end lqg;

  end Design;

encapsulated package Plot "Functions to plot state space system responses"

    encapsulated function polesAndZeros
      "Plot poles (i.e. eigenvalues) and/or invariant zeros of a state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss "Linear system in state space form" annotation(Dialog);
      input Boolean poles=true
        "= true, to plot the poles (i.e. the eigenvalues) of ss"                        annotation(choices(__Dymola_checkBox=true));
      input Boolean zeros=true "= true, to plot the (invariant) zeros of ss " annotation(choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(
         defaultDiagram = Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros());

    protected
      Integer nx=size(ss.A, 1);
      Real eval[nx,2];
      Real invZerosRe[:];
      Real invZerosIm[:];
      Complex invZeros[:];
      Plot.Records.Curve curves[2];
      Integer i;
      Plot.Records.Diagram diagram2;
    algorithm
      // Determine eigen values
      if poles then
        eval := Modelica.Math.Matrices.eigenValues(ss.A);
      end if;

      if zeros then
        invZeros := StateSpace.Analysis.invariantZeros(ss);
        invZerosRe := fill(0,size(invZeros,1));
        invZerosIm := fill(0,size(invZeros,1));
        for i in 1:size(invZeros, 1) loop
          invZerosRe[i] := invZeros[i].re;
          invZerosIm[i] := invZeros[i].im;
        end for;
      end if;

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

           annotation (interactive=true, Documentation(info="<html>
<h4><font style=\"color: #008000; \">Syntax</font></h4>
<blockquote><pre>
StateSpace.Plot.<b>polesAndZeros</b>(ss);
   or
diagram = StateSpace.Plot.<b>polesAndZeros</b>(ss, poles=true, zeros=true, plot=true,
                     defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros</a>(),
                     device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>()); 
</pre></blockquote>

<h4><font style=\"color: #008000; \">Description</font></h4>
<p>
This function plots a pole-zero-map of the poles and transmission zeros of a state space system.
The poles are the eigenvalues of the system matrix (eigenvalues(ss.A)). The Boolean inputs
\"poles\" and \"zeros\" define what to plot. If Boolean input \"plot = true\", the pole-zero-map
is plotted. If false, only the diagram is generated and returned as output argument.
The records \"defaultDiagram\" and \"device\" allow to set various layout options and the
size and location of the diagram on the screen.
</p>

<h4><font style=\"color: #008000; \">Example</font></h4>

<p>
The example <a href=\"Modelica://Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros\">
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
 
<blockquote><img src=\"modelica://Modelica_LinearSystems2/Extras/Images/PolesAndZeros.png\"/> </blockquote>

</html>"));
    end polesAndZeros;

    encapsulated function bodeSISO
      "Plot bode plot of the corresponding transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

      input StateSpace ss "state space system";
      input Integer iu=1 "index of input";
      input Integer iy=1 "index of output";
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min=0.1
        "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max=10
        "Maximum frequency value, if autoRange = false";

      input Boolean magnitude=true "= true, to plot the magnitude of tf" 
                                                                        annotation(choices(__Dymola_checkBox=true));
      input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

    protected
      ZerosAndPoles zp "ZP-Transfer functions to be plotted";
      StateSpace ss_siso(
        redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1),1],
        redeclare Real C[1,size(ss.C, 2)],
        redeclare Real D[1,1]);

    algorithm
      assert(iu <= size(ss.B, 2) and iu > 0, "index for input is " + String(iu) + " which is not in [1, "
         + String(size(ss.B, 2)) + "].");
      assert(iy <= size(ss.C, 1) and iy > 0, "index for output is " + String(iy) + " which is not in [1, "
         + String(size(ss.C, 1)) + "].");
      ss_siso := StateSpace(
        A=ss.A,
        B=matrix(ss.B[:, iu]),
        C=transpose(matrix(ss.C[iy, :])),
        D=matrix(ss.D[iy, iu]));
      zp := StateSpace.Conversion.toZerosAndPoles(ss_siso);

      ZerosAndPoles.Plot.bode(
        zp,
        nPoints,
        autoRange,
        f_min,
        f_max,
        magnitude,
        phase,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
StateSpace.Plot.<b>plotBodeSISO</b>(ss)
   or
StateSpace.Plot.<b>plotBodeSISO</b>(ss, iu, iy, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Plots the bode-diagram of a transfer function.
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>plotBodeSISO</b> plots a bode-diagram of the transfer function corresponding to the behavior of the state space system from iu'th element of the input vector <b>u</b> to the iy'th element of the output vector <b>y</b>.


</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagnitude.png\">
</p>
<p>
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhase.png\">
</p>
<p>


</html> "));
    end bodeSISO;

    encapsulated function bodeMIMO
      "Plot bode plot of all transfer functions, corresponding to the state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss;
      input Integer nPoints(min=2) = 200 "Number of points";
      input Boolean autoRange[:,:]=fill(
          true,
          size(ss.C, 1),
          size(ss.B, 2))
        "= true, if abszissa range is automatically determined";
      input Modelica.SIunits.Frequency f_min[:,:]=fill(
          0.1,
          size(ss.C, 1),
          size(ss.B, 2)) "Minimum frequency value, if autoRange = false";
      input Modelica.SIunits.Frequency f_max[:,:]=fill(
          10,
          size(ss.C, 1),
          size(ss.B, 2)) "Maximum frequency value, if autoRange = false";
      input Boolean magnitude=true "= true, to plot the magnitude of tf" 
                                                                        annotation(choices(__Dymola_checkBox=true));
      input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
            Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

    protected
      TransferFunction tf[size(ss.C, 1),size(ss.B, 2)]
        "Transfer functions to be plotted";
      Plot.Records.Diagram diagram2=defaultDiagram;
      String yNames[size(ss.C, 1)];
      String uNames[size(ss.B, 2)];

    algorithm
     // generate headings
      for i1 in 1:size(ss.B, 2) loop
        uNames[i1] := if ss.uNames[i1] == "" then "u" + String(i1) else ss.uNames[i1];
      end for;
      for i1 in 1:size(ss.C, 1) loop
        yNames[i1] := if ss.yNames[i1] == "" then "y" + String(i1) else ss.yNames[i1];
      end for;

      tf := StateSpace.Conversion.toTransferFunctionMIMO(ss);

      for i1 in 1:size(ss.C, 1) loop
        for i2 in 1:size(ss.B, 2) loop
          diagram2.heading := defaultDiagram.heading +"  "+ uNames[i2] + " -> " + yNames[i1];
          TransferFunction.Plot.bode(
            tf[i1, i2],
            nPoints,
            autoRange[i1, i2],
            f_min[i1, i2],
            f_max[i1, i2],
            magnitude,
            phase,
            defaultDiagram=diagram2,
            device=device);
        end for;
      end for;

      annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
StateSpace.Plot.<b>plotBodeMIMO</b>(ss)
   or
StateSpace.Plot.<b>plotBodeMIMO</b>(ss, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>
</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagnitude.png\">
</p>
<p>
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhase.png\">
</p>
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagnitude2.png\">
</p>
<p>
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhase2.png\">
</p>


</html> "));
    end bodeMIMO;

    encapsulated function timeResponse
      "Plot the time response of the system. The response type is selectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss;
      input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";

      input Modelica_LinearSystems2.Types.TimeResponse response=
          Modelica_LinearSystems2.Types.TimeResponse.Step;

      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                         annotation(Dialog,choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
            heading="Time response"));

    protected
      Plot.Records.Curve curve;
      Integer i1;
      Integer i2;
      Plot.Records.Diagram diagram2[size(ss.C, 1)];

      Real y[:,size(ss.C, 1),if response == TimeResponse.Initial then 1 else size(ss.B,2)]
        "Output response: (number of samples) x (number of outputs) x (number of inputs)";
      Real t[:] "Time vector: (number of samples)";
      Real x[:,size(ss.A, 1),if response == TimeResponse.Initial then 1 else size(ss.B,2)]
        "State trajectories: (number of samples) x (number of states) x (number of inputs)";
      String yNames[size(ss.C, 1)];
      String uNames[size(ss.B, 2)];
      Integer loops=if response == TimeResponse.Initial then 1 else size(ss.B,2);

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
          diagram2[i1].heading := if response == TimeResponse.Initial then defaultDiagram.heading +" "+ yNames[i1] else defaultDiagram.heading + "  " + uNames[i2] + " -> " + yNames[i1];
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

      annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>timeResponse</b>(ss);
or
StateSpace.Plot.<b>timeResponse</b>(ss, dt, tSpan,response, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>


<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>timeResponse</b> plots the time response of a state space system. The character of the time response if defined by the input <tt>response</tt>, i.e. Impulse, Step, Ramp, or Initial. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.plotImpulse\">plotImpulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/timeResponseSS.png\">
</p>
<p>
</p>


</html> "));
    end timeResponse;

    encapsulated function impulse "Impulse response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss;
      input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                         annotation(Dialog,choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response"));

    protected
      input Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Impulse
        "type of time response";
    algorithm

      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
        ss=ss,
        dt=dt,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>impulse</b>(ss);
or
StateSpace.Plot.<b>impulse</b>(ss, dt, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>plotImpulse</b> plots the impulse responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/impulseResponseSS.png\">
</p>
<p>
</p>


</html> "));
    end impulse;

    encapsulated function step "Step response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss;
      input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                         annotation(Dialog,choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response"));

      input Modelica_LinearSystems2.Types.TimeResponse response=
          Modelica_LinearSystems2.Types.TimeResponse.Step
        "type of time response";

    algorithm
      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
        ss=ss,
        dt=dt,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>step</b>(ss);
or
StateSpace.Plot.<b>step</b>(ss, dt, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>step</b> plots the step responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.initial\">initial</a>.





</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/stepResponseSS.png\">
</p>
<p>
</p>


</html> "));
    end step;

    encapsulated function ramp "Ramp response plot"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input StateSpace ss;
    input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";
    input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

    input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                         annotation(Dialog,choices(__Dymola_checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
          heading="Ramp response"));

    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Ramp "type of time response";

    algorithm
    Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
          ss=ss,
          dt=dt,
          tSpan=tSpan,
          response=response,
          x0=x0,
          defaultDiagram=defaultDiagram,
          device=device);

    annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>ramp</b>(ss);
or
StateSpace.Plot.<b>ramp</b>(ss, dt, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>ramp</b> plots the ramp responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/rampResponseSS.png\">
</p>
<p>
</p>


</html> "));
    end ramp;

    encapsulated function initialResponse "Initial condition response plot"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

      input StateSpace ss;
      input Real dt=0 "Sample time [s]";
      input Real tSpan=0 "Simulation time span [s]";
      input Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";

      input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                         annotation(Dialog,choices(__Dymola_checkBox=true));

      extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
            heading="Initial response"));

      input Modelica_LinearSystems2.Types.TimeResponse response=
          Modelica_LinearSystems2.Types.TimeResponse.Initial
        "type of time response";
    algorithm

      Modelica_LinearSystems2.StateSpace.Plot.timeResponse(
        ss=ss,
        dt=dt,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
StateSpace.Plot.<b>initial</b>(ss);
or
StateSpace.Plot.<b>initial</b>(ss, dt, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initial</b> plots the initial responses of a state space system for the initial state vector x0 for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.step\">step</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Plot.ramp\">ramp</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
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

</p>
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/initialResponseSS.png\">
</p>
<p>
</p>


</html> "));
    end initialResponse;

end Plot;

encapsulated package Conversion
    "Conversion functions from StateSpace into TransferFunction or ZerosAndPoles representations"

  encapsulated function toZerosAndPoles
      "Generate a zeros-and-poles representation from a SISO state space representation"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.Math.Complex;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss "StateSpace object";
  //protected
  //  input Boolean cancel=true "false to hinder cancellation";// is not fully realized
    public
    output ZerosAndPoles zp;

    protected
    StateSpace ssm= if size(ss.A,1)>0 then StateSpace.Transformation.toIrreducibleForm(ss) else StateSpace(ss.D[1,1]);
    Complex poles[:];
    Complex zeros[:];

    Real gain;

    Complex frequency;
    Complex Gs;
    Real As[:,:];
    Real pk;
    Integer i;
    Integer k;
    Boolean h;
  Real v;

  algorithm
    if Modelica.Math.Vectors.length(ssm.B[:, 1]) > 0 and 
        Modelica.Math.Vectors.length(ssm.C[1, :]) > 0 then

      poles := Complex.Internal.eigenValues_dhseqr(ssm.A);//ssm.A is of upper Hessenberg form

     zeros := StateSpace.Analysis.invariantZeros(ssm);

      if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
        assert(size(ss.B, 2) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
           + String(size(ss.B, 2)) + " instead of 1");
        assert(size(ss.C, 1) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
           + String(size(ss.C, 1)) + " instead of 1");
      end if;
      zp := ZerosAndPoles(
          z=zeros,
          p=poles,
          k=1);

      v := sum(cat(1, zeros[:].re,  poles[:].re))/max(size(zeros,1)+size(poles,1),1)+13/19;
      frequency := Complex(v)*19/17;

      Gs := ZerosAndPoles.Analysis.evaluate(zp, frequency);

      As := -ssm.A;
      for i in 1:size(As, 1) loop
        As[i, i] := As[i, i] + frequency.re;
      end for;

      pk := StateSpace.Internal.partialGain(As, ssm.B[:, 1]);
      gain := (ssm.C[1, size(As, 1)]*pk + ss.D[1, 1])/Gs.re;

      zp := ZerosAndPoles(
          z=zeros,
          p=poles,
          k=gain);

    else
      zp := ZerosAndPoles(
          z=fill(Complex(0), 0),
          p=fill(Complex(0), 0),
          k=scalar(ss.D));

    end if;
    zp.uName := ss.uNames[1];
    zp.yName := ss.yNames[1];

    annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> StateSpace.Conversion.<b>toZerosAndPoles</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a ZerosAndPoles record
 <blockquote><pre>
                 product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
        zp = k*---------------------------------------------------------
                product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>of a system from state space representation using the transformation algorithm described in [1].
<br>
The uncontrollable and unobservable parts are isolated and the eigenvalues and invariant zeros of the controllable and observable sub system are calculated.


<h4><font color=\"#008000\">Example</font></h4>
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


<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V.  </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control,
vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>

</html> "));
  end toZerosAndPoles;

  function toTransferFunction
      "Generate a transfer function from a SISO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

    input StateSpace ss "StateSpace object";

    output TransferFunction tf;

    protected
    ZerosAndPoles zp;

  algorithm
    zp := toZerosAndPoles(ss);
    tf := ZerosAndPoles.Conversion.toTransferFunction(zp);

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> StateSpace.Conversion.<b>toTransferFunction</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a TransferFunction record
<blockquote><pre>
           n(s)     b0 + b1*s + ... + bn*s^n
   tf = -------- = -------------------------- 
           d(s)     a0 + a1*s + ... + an*s^n
 </pre></blockquote>

The algorithm uses <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPoles\">toZerosAndPoles</a> to convert the state space system into a zeros and poles representation first and after that href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction\">ZerosAndPoles.Conversion.toTransferFunction</a> to generate the transfer function.



<h4><font color=\"#008000\">Example</font></h4>
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




</html> "));
  end toTransferFunction;

encapsulated function toZerosAndPolesMIMO
      "Generate a zeros-and-poles representation from a MIMO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

  input StateSpace ss "StateSpace object";

  output ZerosAndPoles zp[size(ss.C, 1),size(ss.B, 2)];

    protected
  StateSpace ss_siso(
    redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
    redeclare Real B[size(ss.B, 1),1],
    redeclare Real C[1,size(ss.C, 2)],
    redeclare Real D[1,1]);

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
          zp[ic, ib] := StateSpace.Conversion.toZerosAndPoles(ss_siso);
     end for;
  end for;
  annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> StateSpace.Conversion.<b>toZerosAndPolesMIMO</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a matrix of ZerosAndPoles records
 <blockquote><pre>
                 product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
        zp = k*---------------------------------------------------------
                product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
of a system from state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each input-output pair.


<h4><font color=\"#008000\">Example</font></h4>
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
i.e.
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



<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V.  </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control,
vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>

</html> "));
end toZerosAndPolesMIMO;

function toTransferFunctionMIMO
      "Generate a transfer function of a MIMO system from state space representation"
      import Modelica_LinearSystems2;

      import Modelica;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;

  input StateSpace ss "StateSpace object";

  output TransferFunction tf[size(ss.C, 1),size(ss.B, 2)]
        "Matrix of transfer function objects";

    protected
  ZerosAndPoles zp[:,:];
  parameter Integer m=size(ss.B, 2);
  parameter Integer p=size(ss.C, 1);

algorithm
  zp := Modelica_LinearSystems2.StateSpace.Conversion.toZerosAndPolesMIMO(ss);
  for i1 in 1:m loop
    for i2 in 1:p loop
      tf[i2, i1] := ZerosAndPoles.Conversion.toTransferFunction(zp[i2, i1]);
    end for;
  end for;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  zp </td><td align=center> =  </td>  <td> StateSpace.Conversion.<b>toTransferFunctionMIMO</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a matrix of TransferFunction records
<blockquote><pre>
           n(s)     b0 + b1*s + ... + bn*s^n
   tf = -------- = -------------------------- 
           d(s)     a0 + a1*s + ... + an*s^n
 </pre></blockquote>
with repetitive application of <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Conversion.toTransferFunction\">Conversion.toTransferFunction</a>


<h4><font color=\"#008000\">Example</font></h4>
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
i.e.
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




</html> "));
end toTransferFunctionMIMO;

end Conversion;

  encapsulated package Transformation "State Space similarity transformations"

      encapsulated function toSimilarForm
      "Perform the similarity transformation z = Tx (or x = inv(T)z) which leads to Az=T*A*inv(T), Bz=T*B, Cz=C*inv(T), Dz=D (or Az=inv(T)*A*T, Bz=inv(T)B, Cz=C*T, Dz=D)"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input StateSpace ss "state space system";
        input Real T[size(ss.A, 2),size(ss.A, 1)] = identity(size(ss.A,1))
        "Transformation matrix";
        input Boolean inverted=false
        "false (default) for transformation z = Tx, true for x = Tz"   annotation(choices(checkBox=true));

        output StateSpace tss(
          redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
          redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

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
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tss </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toSimilarForm</b>(ss, T, inverted)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>toSimilarForm</b> calculates a similar state space system, i.e.
<blockquote><pre>
   der(z) = T*A*inv(T)*z + T*B*u
        y = C*T*z + D*u
</pre></blockquote> 
if inverted==false and 
<blockquote><pre>
   der(z) = inv(T)*A*T*z + inv(T)*B*u
        y = C*inv(T)*z + D*u
</pre></blockquote> 
if inverted=true. Matrix T has to be invertible. The transformed system has the same eigenvalues. See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.analysis\">analysis</a>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
      end toSimilarForm;

      encapsulated function toObservabilityForm
      "Perform the similarity transformation to the obervabillity canonical form"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input StateSpace ss "state space system";
        output StateSpace tss(
          redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
          redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

    protected
        Real V[size(ss.A, 1),size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

        Integer nx=size(ss.A, 1);

      algorithm
        assert(size(ss.C, 1) == 1 and size(ss.B, 2) == 1,
          "Calculation of controllable form fails for systems with more than 1 inputs or outputs");
        assert(Modelica_LinearSystems2.StateSpace.Analysis.isObservable(ss),
          "transformation ist not realizable since the system ist not obersvable");

        V[:, 1] := Modelica.Math.Matrices.solve(StateSpace.Analysis.observabilityMatrix(ss), vector([
          fill(
            0,
            1,
            nx - 1),1]));

        for i in 2:nx loop
          V[:, i] := ss.A*V[:, i - 1];
        end for;

        tss := StateSpace.Transformation.toSimilarForm(ss, V,inverted=true);

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tss </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toObservabilityForm</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>toObservabilityForm</b> computes the observability form of a SISO state space system, i.e.
<blockquote><pre>
   tss:
   der(z) = inv(T)*A*T*z + inv(T)*B*u
        y = C*inv(T)*z + D*u
</pre></blockquote>
with
<blockquote><pre>
   T = [C; C*A; ...; C*A^(n-1)]
</pre></blockquote>
is the observability matrix of the original state space system.
In comparison to the corresponding transfer function
<blockquote><pre>
           b0 + b1*s + ... + bn*s^n
   G(s) = --------------------------
           a0 + a1*s + ... + an*s^n
 </pre></blockquote>
the canonical observability form is
<blockquote><pre>
       | 0   0   ...   0   -a0   |        | b0   - a0*bn   |
       | 1   0   ...   0   -a1   |        | b1   - a1*bn   |
   A = | 0   1   ...   0   -a2   |,   B = |     ...        |
       |... ...  ...  ...  -a3   |        | bn-2 - an-2*bn |
       | 0  ...  ...   1   -an-1 |        | bn-1 - an-1*bn |
       
   C = [0, 0, ..., 1],                D = [bn]
 </pre></blockquote>




Matrix T has to be invertible, i.e. the system has to be observable. The transformed system has the same eigenvalues. See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>, <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toControllabilityForm\">toControllabilityForm</a>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
      end toObservabilityForm;

      encapsulated function toControllabilityForm
      "Perform the similarity transformation to the controllability canonical form"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Streams;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input StateSpace ss "state space system";
        output StateSpace tss(
          redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
          redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

    protected
        Real V[size(ss.A, 1),size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

        Integer nx=size(ss.A, 1);

      algorithm
        assert(size(ss.C, 1) == 1 and size(ss.B, 2) == 1,
          "Calculation of controllable form fails for systems with more than 1 inputs or outputs");
        assert(Modelica_LinearSystems2.StateSpace.Analysis.isControllable(ss),
          "transformation ist not realizable since the system ist not controllable");

        V[1, :] := Modelica.Math.Matrices.solve(transpose(StateSpace.Analysis.controllabilityMatrix(ss)),
          vector([fill(
            0,
            nx - 1,
            1); 1]));

        for i in 2:nx loop
          V[i, :] := V[i - 1, :]*ss.A;
        end for;

        tss := StateSpace.Transformation.toSimilarForm(
                                ss, V);

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tss </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toControllabilityForm</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>toControllabilityForm</b> computes the controllability form of a SISO state space system, i.e.
<blockquote><pre>
   tss:
   der(z) = T*A*inv(T)*z + T*B*u
        y = C*T*z + D*u
</pre></blockquote>
with
<blockquote><pre>
   T = [B, A*B,..., A^(n-1)*B]
</pre></blockquote>
is the observability matrix of the original state space system.
In comparison to the corresponding transfer function
<blockquote><pre>
           b0 + b1*s + ... + bn*s^n
   G(s) = --------------------------
           a0 + a1*s + ... + an*s^n
 </pre></blockquote>
the canonical observability form is
<blockquote><pre>
       | 0   1   0   ...   0     0   |                        | 0 |
       |  0     0     0    0     0   |                        | 0 |
   A = | ...   ...   ...  ...   ...  |,                   B = |...|
       |  0     0     0    0     0   |                        | 0 |
       | -a0   -a1   -a2  ...  -an-1 |                        | 1 |
       
   C = [ b0 - bn*a0, b1 - bn*a1, ..., bn-1 - bn*an-1],    D = [bn]
 </pre></blockquote>




Matrix T has to be invertible, i.e. the system has to be controllable. The transformed system has the same eigenvalues. See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>, <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toObservabilityForm\">toObservabilityForm</a>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
      end toControllabilityForm;

      encapsulated function toDiagonalForm
      "Perform the similarity transformation with the (real) inverse right eigenvector matrix of the system, that lead to the Jordan canonical form for single eigenvalues"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

        input StateSpace ss "state space system";
        output StateSpace tss(
          redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
          redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

    protected
        Real V[size(ss.A, 1),size(ss.A, 1)]
        "Matrix of the right eigenvectors of the matrix ss.A";

      algorithm
        (,,,V,) := LAPACK.dgeev(ss.A);

        tss := StateSpace.Transformation.toSimilarForm(ss, V, inverted=true);

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tss </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toDiagonalForm</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>toDiagonalForm</b> computes the diagonal form of a SISO state space system, i.e.
<blockquote><pre>
   tss:
   der(z) = inv(T)*A*T*z + inv(T)*B*u
        y = C*inv(T)*z + D*u
</pre></blockquote>


Matrix T has to be diagonalizable, i.e. the algebraic and geometric multiplicities of an eigenvalue must coincide. The diagonal entries of the new system matrix tss.<b>A</b> are the eigenvalues off the systemmatrix ss.<b>A</b>. See also <a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Transformation.toSimilarForm\">toSimilarForm</a>.

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
      end toDiagonalForm;

      encapsulated function toIrreducibleForm
      "Calculate a minimal controllable and observable block Hessenberg realization of a given SISO state-space representation "

       // test of SISO has to be added
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica;

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
        Integer rankQ2=ssm2.r;
    public
        output StateSpace ssm3(
          redeclare Real A[rankQ2,rankQ2],
          redeclare Real B[rankQ2,size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),rankQ2],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);
      algorithm
        ssm3 := StateSpace(
            A=transpose(ssm2.A[nx2 - rankQ2 + 1:nx2, nx2 - rankQ2 + 1:nx2]),
            B=transpose(ssm2.C[:, nx2 - rankQ2 + 1:nx2]),
            C=transpose(ssm2.B[nx2 - rankQ2 + 1:nx2, :]),
            D=(ssm2.D));

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  tss </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toIrreducibleForm</b>(ss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function calculates a minimal controllable and observable block Hessenberg realization for a given state-space representation.
Therefore, all uncontrollable and unobservable modes are removed by performing orthogonal similarity transformations as described in [1].
<p>
This function is called to compute transfer functions of state space representations as described in [1]. Look at [1] for further details
<h4><font color=\"#008000\">Example</font></h4>
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
      D=[0]
)
</pre></blockquote>
<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V. </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control, vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>
</html> "));
      end toIrreducibleForm;

      encapsulated function toStaircaseForm
      "Transforms a state space system to upper staircase form"
        import Modelica;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.StateSpace;

        input StateSpace ss "state space system";
        input Modelica_LinearSystems2.Types.StaircaseMethod method=
            Modelica_LinearSystems2.Types.StaircaseMethod.SVD;

        output StateSpace ss_sc(
          redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
          redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
          redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
          redeclare Real D[size(ss.D, 1),size(ss.D, 2)]);

      algorithm
      if method == Modelica_LinearSystems2.Types.StaircaseMethod.SVD then
        (,ss_sc) := StateSpace.Internal.staircaseSVD(ss);
      else
        (,ss_sc) := Modelica_LinearSystems2.StateSpace.Internal.staircaseQR(ss);
        end if;

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss_sc </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>toStaircaseForm</b>(ss, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>toStaircaseForm</b> computes the upper staircase form state space system.


<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
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
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(inputIndex, 1)],
      redeclare Real C[size(outputIndex, 1),size(ss.C, 2)],
      redeclare Real D[size(outputIndex, 1),size(inputIndex, 1)])
        "Subsystem state space record";

  algorithm
    subSc.A := ss.A;
    subSc.B := ss.B[:, inputIndex];
    subSc.C := ss.C[outputIndex, :];
    subSc.D := ss.D[outputIndex, inputIndex];

        annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  subsystem </td><td align=center> =  </td>  <td> StateSpace.Transformation.<b>extract</b>(ss, outputIndex, inputIndex)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>extract</b> computes the subsystem of a state space system corresponding to the indices in outputIndex and inputIndex, i.e.
<blockquote><pre>
  subsystem.A = ss.A;
  subsystem.B = ss.B[:, inputIndex];
  subsystem.C = ss.C[outputIndex, :];
  subsystem.D = ss.D[outputIndex, inputIndex];</pre></blockquote>

<h4><font color=\"#008000\">Example</font></h4>
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


</html> "));
  end extract;

  end Transformation;

encapsulated package Import
    "Utilitiy functions to import StaeSpace representations"

  encapsulated function fromFile "Read a StateSpace data record from mat-file"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace.Internal;
    input String fileName="dslin.mat"
        "Name of the state space system data file"     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    input String matrixName="ABCD" "Name of the state space system matrix"    annotation(Dialog);
    protected
    input Integer xuy[3]=Internal.readSystemDimension(fileName, matrixName);
    input Integer nx=xuy[1];
    input Integer nu=xuy[2];
    input Integer ny=xuy[3];

    public
    output StateSpace result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "= model read from file";

    protected
    Real ABCD[nx + ny,nx + nu]=Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
          fileName,
          matrixName,
          nx + ny,
          nx + nu);

  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    Modelica.Utilities.Streams.print("StateSpace record loaded from file: \""
       + fileName + "\"");

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> StateSpace.Import.<b>fromFile</b>(fileName, matrixName)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Reads and loads a state space system from a mat-file <tt>fileName</tt>. The file must contain the matrix [A, B; C, D] named matrixName and the integer nx representing the order of the system, i.e. the number of rows of the square matrix A.

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
     

<b>algorithm</b>
  ss:=Modelica_LinearSystems2.StateSpace.Import.fromFile(\"stateSpace.mat\", \"ABCD\");
//  ss=StateSpace(
      A=[-1, 0, 0; 0, -2, 0; 0, 0, -3],
      B=[1; 1; 0],
      C=[1, 1, 1],
      D=[0])


</pre></blockquote>


</html> "));
  end fromFile;

  function fromModel
      "Generate a StateSpace data record by linearization of a model"

      import Modelica;
      import Modelica_LinearSystems2.StateSpace;

    input String modelName "Name of the Modelica model" annotation(Dialog(translatedModel));
    input Real T_linearize=0
        "point in time of simulation to linearize the model";
    input String fileName="dslin" "Name of the result file";
    protected
    String fileName2=fileName + ".mat";
    Boolean OK1 = simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
    Boolean OK2 = importInitial("dsfinal.txt");
    Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize+1);

    Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
    Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
    Integer nx=integer(nxMat[1, 1]);
    Integer nu=ABCDsizes[2] - nx;
    Integer ny=ABCDsizes[1] - nx;
    Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
    String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);
    public
    output StateSpace result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "= model linearized at initial point";

  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    result.uNames := xuyName[nx + 1:nx + nu];
    result.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    result.xNames := xuyName[1:nx];

          annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> StateSpace.Import.<b>fromModel</b>(modelName, T_linearize, fileName)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Generate a StateSpace data record by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation. The result of linearization is transformed into a StateSpace record.

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   String modelName = \"Modelica_LinearSystems2.Examples.Utilities.DoublePendulum\"; 
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



</html> 
"));
  end fromModel;

    annotation (Documentation(info="<html>
</html>"));
end Import;

encapsulated package Internal
    "Internal library of record StateSpace (should not be directly used by user)"
    import Modelica;
    import Modelica_LinearSystems2;
  extends Modelica.Icons.Library;

  encapsulated function isSISO
      "To check a state space system to be SISO (or not)"

      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

    output Boolean isSISO;
  algorithm
    isSISO := size(ss.B, 2) == 1 and size(ss.C, 1) == 1;
    annotation (Documentation(info="<html>
 
 
</html>"));
  end isSISO;

encapsulated function invariantZeros2
      "Compute invariant zeros of linear SISO state space system with a generalized system matrix [A, B, C, D] which is of upper Hessenberg form"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;
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
  Real A[nx + ny,nx + nu]=[ss.A,ss.B; ss.C,ss.D];
  Real B[nx + ny,nx + nu]=[identity(nx),zeros(nx, nu); zeros(ny, nx + nu)];
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
inputs (= "
          + String(nu) + ") = number of outputs (= " + String(ny) + ")
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
      z[j].re := if abs(alphaReal[i]) >= normB*1e-12 then alphaReal[i]/beta[i] else 
              0;
      z[j].im := if abs(alphaImag[i]) >= normB*1e-12 then alphaImag[i]/beta[i] else 
              0;
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
The advantage of this function in comparision to the general invariantZeros function
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

function characterizeEigenvalue
      "Check stability, stabilizability, controllability, observability nad detectability of the single poles"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.Eigenvalue;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Internal;

  input StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=fill(
        0,
        0,
        0),
      B=fill(
        0,
        0,
        0),
      C=fill(
        0,
        0,
        0),
      D=fill(
        0,
        0,
        0));
  input Eigenvalue evin[:];
  output Eigenvalue ev[size(ss.A, 1)];

    protected
  Real cPoles[:,2] "controllable poles";
  Real ncPoles[:,2] "uncontrollable poles";
  Real poles[size(ss.A, 1),2] "controllable and uncontrollable poles";
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
  Complex j = Modelica_LinearSystems2.Math.Complex.j();

algorithm
  for i in 1:nx loop
    ev[i].ev := evin[i].ev;
  end for;

  (cPoles,ncPoles,poles) := StateSpace.Internal.controllablePoles(ss);
  n_c := size(cPoles, 1);
  for i in 1:n_c loop
    cdummy[i].ev := cPoles[i, 1]+cPoles[i, 2]*j;
    cdummy[i].isControllable := true;
  end for;
  for i in 1:size(ncPoles, 1) loop
    cdummy[n_c + i].ev := ncPoles[i, 1]+ncPoles[i, 2]*j;
    cdummy[n_c + i].isControllable := false;
  end for;

// controllable poles of the trnasposed system are the observable poles
  (cPoles,ncPoles,poles) := StateSpace.Internal.controllablePoles(sst);
  n_c := size(cPoles, 1);
  for i in 1:n_c loop
    odummy[i].ev := cPoles[i, 1]+cPoles[i, 2]*j;
    odummy[i].isObservable := true;
  end for;
  for i in 1:size(ncPoles, 1) loop
    odummy[n_c + i].ev := ncPoles[i, 1]+ncPoles[i, 2]*j;
    odummy[n_c + i].isObservable := false;
  end for;

// using ev.imag as an flag to mark the eigenvalues
  for i in 1:nx loop
    absVector := fill(1e50, nx);
    for ii in 1:nx loop
      if not odummy[ii].imag then
        absVector[ii] := Complex.'abs'(ev[i].ev - odummy[ii].ev);
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
        absVector[ii] := Complex.'abs'(ev[i].ev - cdummy[ii].ev);
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

  encapsulated function isStabilizableSISO
      "To check wether a SISO system is stabliziable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    protected
    input Modelica_LinearSystems2.Internal.StateSpaceR ssm=
        StateSpace.Internal.cntrHessenberg(ss);
    public
    output Boolean stabilizable;

    protected
    Complex evd[:]=fill(Complex(0), size(ss.A, 1) - ssm.r);

  algorithm
    if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
      assert(size(ss.B, 2) == 1,
        "A SISO-system is expected as input\n but the number of inputs is " +
        String(size(ss.B, 2)) + " instead of 1");
      assert(size(ss.C, 1) == 1,
        " A SISO-system is expected as input\n but the number of outputs is " +
        String(size(ss.C, 1)) + " instead of 1");
    end if;
    evd := Complex.eigenValues(ssm.A[ssm.r + 1:size(ss.A, 1), ssm.r + 1:size(ss.A,
      1)]);
    stabilizable := true;

    if size(ss.A, 1) > ssm.r then
       for i1 in 1:size(evd, 1) loop
        stabilizable := stabilizable and evd[i1].re < 0;
      end for;
    end if;

      annotation (Documentation(info="<html>
This function checks whether a SISO state space system is stabilizable or not.
<p>
A system is stabilizable for the continuous-time case if all of the uncontrollable eigenvalues have neagtive real part
or for the discrete-time case if all of the uncontrollable eigenvalues are in the complex unit circle respectively.
Hence, a controllable system is always stabilizable of course.
<p>
To check stabilizability, ths system is transformed to to upper controller Hessenberg form
<blockquote><pre>
               | *   *   ...   ...    * |               | * |
               | *   *   ...   ...    * |               | 0 |
 <b>Q</b>*<b>A</b>*<b>Q</b> ' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>b</b> = <b>q</b> = | . |,   <b>c</b>*<b>Q</b> = ( *, ..., * )
               | .   .    .     .     . |               | . |
               | 0  ...   0     *     * |               | 0 |
 
</pre>
</blockquote>
The system can be partitioned to 
 
<blockquote><pre>
<b>H</b>=[<b>H</b>11,<b>H</b>12; <b>H</b>21, <b>H</b>22], <b>q</b>=[<b>q</b>1;<b>0</b>],
</pre>
</blockquote
where the pair (<b>H</b>11, <b>q</b>1) contains the controllable part of the system, that is, rank(<b>H</b>) = rank(<b>H</b>11). For
stabilizability the <b>H</b>22 has to be stable.
 
 
 
 
</html>"));
  end isStabilizableSISO;

  encapsulated function isStabilizableMIMO
      "To check wether a MIMO system is stabliziable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;

      input StateSpace ss;

      output Boolean stabilizable;

    protected
      Complex evnc[:] "complex vector of uncontrollable poles";
      Real cPoles[:,2] "controllable poles";
      Real ncPoles[:,2] "uncontrollable poles";
      Real poles[size(ss.A, 1),2] "controllable and uncontrollable poles";
      Complex j=Modelica_LinearSystems2.Math.Complex.j();

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

  encapsulated function isDetectableSISO
      "To check wether a SISO system is detectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    protected
    input Modelica_LinearSystems2.Internal.StateSpaceR ssm=
        StateSpace.Internal.cntrHessenberg(
        StateSpace.Internal.transposeStateSpace(ss));
    public
    output Boolean detectable;
    protected
    Complex evd[:]=fill(Complex(0), size(ss.A, 1) - ssm.r);

  algorithm
    if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
      assert(size(ss.B, 2) == 1,
        "A SISO-system is expected as input\n but the number of inputs is " +
        String(size(ss.B, 2)) + " instead of 1");
      assert(size(ss.C, 1) == 1,
        " A SISO-system is expected as input\n but the number of outputs is " +
        String(size(ss.C, 1)) + " instead of 1");
    end if;
    evd := Complex.eigenValues(ssm.A[ssm.r + 1:size(ss.A, 1), ssm.r + 1:size(ss.A,1)]);

    if size(ss.A, 1) == ssm.r then
      detectable := true;

    else
      detectable := true;
      for i1 in 1:size(evd, 1) loop
        detectable := detectable and evd[i1].re < 0;
      end for;
    end if;

      annotation (Documentation(info="<html>
This function checks whether a SISO state space system is detectable or not.
<p>
A system is detectable for the continuous-time case if all of the unobservable eigenvalues have neagtive real part
or for the discrete-time case if all of the unobservable eigenvalues are in the complex unit circle respectively.
Hence, a oberservable system is always detectable of course.
<p>
As observability is a dual concept of controllability, the concept of detectability is dual to stabilizability, that is,
a system is detectable if the pair (<b>A</b>', <b>C</b>') is stabilizable. Therefore, the same algorithm to check stabilizability
are applied to the dual pair (<b>A</b>', <b>C</b>') of the system:
<p>
To check stabilizability (see Modelica_LinearSystems2.StateSpace.Analysis.isStabilizable) , ths system is transformed to to upper controller Hessenberg form
<blockquote><pre>
               | *   *   ...   ...    * |               | * |
               | *   *   ...   ...    * |               | 0 |
 <b>Q</b>*<b>A</b>*<b>Q</b> ' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>b</b> = <b>q</b> = | . |,   <b>c</b>*<b>Q</b> = ( *, ..., * )
               | .   .    .     .     . |               | . |
               | 0  ...   0     *     * |               | 0 |
 
</pre>
</blockquote>
The system can be partitioned to 
 
<blockquote><pre>
<b>H</b>=[<b>H</b>11,<b>H</b>12; <b>H</b>21, <b>H</b>22], <b>q</b>=[<b>q</b>1;<b>0</b>],
</pre>
</blockquote
where the pair (<b>H</b>11, <b>q</b>1) contains the controllable part of the system, that is, rank(<b>H</b>) = rank(<b>H</b>11). For
stabilizability the <b>H</b>22 has to be stable.
</p>
 
</html>"));
  end isDetectableSISO;

  encapsulated function isDetectableMIMO
      "To check wether a MIMO system is detectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;

      input StateSpace ss;

      output Boolean detectable;

    protected
      StateSpace sst=StateSpace.Internal.transposeStateSpace(ss);

      Complex evnd[:] "complex vector of uncontrollable poles";
      Real dPoles[:,2] "controllable poles";
      Real ndPoles[:,2] "uncontrollable poles";
      Real poles[size(ss.A, 1),2] "controllable and uncontrollable poles";
      Complex j=Modelica_LinearSystems2.Math.Complex.j();

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
 
 
 
 
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics));
  end isDetectableMIMO;

  encapsulated function isObservableSISO
      "To check wether a SISO system is observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

     input StateSpace ss;

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
This function is to calculate whether a SISO state space system is observable or not. Therefore, the dual System (A', c', b', d')
it is transformed to upper observer Hessenberg form
<blockquote><pre>
               | *   *   ...   ...    * |             | * |
               | *   *   ...   ...    * |             | . |
 <b>Q</b>*<b>A'</b>*<b>Q</b>' = <b>H</b> = | 0   *   ...   ...    * |,    <b>Q</b>*<b>c'</b> =  | . |,   <b>b'</b>*<b>Q</b> = <b>q</b> = ( 0, ..., 0, * )
               | .   .    .     .     . |             | * |
               | 0  ...   0     *     * |             | * |
  
</pre>
</blockquote>
Note, that
<blockquote><pre>
                         n-1                          n-1
rank(<b>c'</b>; <b>c'*<b>A'</b>; ...; <b>c'</b>*A'</b>   ) = rank(<b>q</b>; <b>q</b>*<b>H</b>; ...; <b>q</b>*<b>H</b>   )
</pre>
</blockquote>
and that
<blockquote><pre>
                  n-1
 (<b>q</b>; <b>H</b>*<b>q</b>; ...; <b>q</b>*<b>H</b>  )
</pre>
</blockquote>
is a lower triangular matrix and has full rank if and only if none of the elements in the diagonal is zero. That is, that neither qn or hi,i-1,   i = 2,..., n   may be zero.
 
 
</html>"));
  end isObservableSISO;

  encapsulated function isControllableSISO
      "To check a SISO system wether it is controllable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

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
      "To check a MIMO system wether it is controllable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      input Modelica_LinearSystems2.Types.StaircaseMethod method=
          Modelica_LinearSystems2.Types.StaircaseMethod.SVD;

      output Boolean controllable;
  algorithm
      assert(method == Modelica_LinearSystems2.Types.StaircaseMethod.SVD or 
        method == Modelica_LinearSystems2.Types.StaircaseMethod.QR, "\nMethods for staircase algorithm are QR factorization or singular value decomposition. Therefore, 
the variable \"method\" in \"Modelica_LinearSystems2.StateSpace.Internal.isControllableMIMO\" has to be qr or svd but is method = "
         + String(method));
      if min(size(ss.B)) == 0 then
        controllable := false;
      else
        if method == Modelica_LinearSystems2.Types.StaircaseMethod.QR then
          controllable := StateSpace.Internal.staircaseQR(ss);
        else
          controllable := StateSpace.Internal.staircaseSVD(ss);
        end if;
      end if;
      annotation (Documentation(info="<html>
 
 
 
</html>"));
  end isControllableMIMO;

  encapsulated function isObservableMIMO
      "To check a MIMO system wether it is observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

      input StateSpace ss;
      input Modelica_LinearSystems2.Types.StaircaseMethod method=
          Modelica_LinearSystems2.Types.StaircaseMethod.SVD;

    protected
      StateSpace ss2=StateSpace.Internal.transposeStateSpace(ss);

    public
      output Boolean observable;
  algorithm
      assert(method == Modelica_LinearSystems2.Types.StaircaseMethod.SVD or 
        method == Modelica_LinearSystems2.Types.StaircaseMethod.QR, "\nMethods for staircase algorithm are QR factorization or singular value decomposition. Therefore, 
the variable \"method\" in \"Modelica_LinearSystems2.StateSpace.Internal.isControllableMIMO\" has to be qr or svd but is method = "
         + String(method));
      if min(size(ss.C)) == 0 then
        observable := false;
      else
        if method == Modelica_LinearSystems2.Types.StaircaseMethod.QR then
          observable := StateSpace.Internal.staircaseQR(ss2);
        else
          observable := StateSpace.Internal.staircaseSVD(ss2);
        end if;

      end if;

      annotation (Documentation(info="<html>
 
 
 
</html>"));
  end isObservableMIMO;

 encapsulated function isControllableAndObservableSISO
      "To check whether a SISO system is controllable and observable"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

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
        "A SISO-system is expected as input\n but the number of inputs is " +
        String(size(ss.B, 2)) + " instead of 1");
      assert(size(ss.C, 1) == 1,
        " A SISO-system is expected as input\n but the number of outputs is "
         + String(size(ss.C, 1)) + " instead of 1");
    end if;
    controllableAndObservable := size(ss.A, 1) == ssm2.r;
 equation

 end isControllableAndObservableSISO;

  encapsulated function readLength_nx
      "Read the order nx of a state space system from a file"

    input String fileName="ss_siso.mat"
        "Name of the state space system data file"     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    output Integer nx;
    protected
    Real nxMat[1,1]=readMatrix(
            fileName,
            "nx",
            1,
            1);
  algorithm
    nx := integer(nxMat[1, 1]);
  end readLength_nx;

  encapsulated function staircaseQR
      "Staircase algorithm to put a state space system to controller Hessenberg form"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

    input StateSpace ss;

    output Boolean isControllable;
    output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
      redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "controllable state space system";
    output Real PP[:,:];
    protected
    Real A[:,:];
    Real B[size(ss.B, 1),size(ss.B, 2)];
    Real C[:,:];
    Real Q[:,:];
    Real Q2[:,:];
    Real R[:,:];

    Real P[:,:];
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

      (Q,R,tau,Q2) := Matrices.QR( ss.B);
      B := [R; zeros(nx - nu, nu)];  // should be the same as transopse(Q2)*ss.B

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

        (Q,R,tau,Q2) := Matrices.QR( A[stairStep + 1:nx, stairStep - rankR +
          1:stairStep]);
        P := [identity(nx - n1),zeros(nx - n1, n1); zeros(n1, nx - n1),Q2];
        PP := transpose(P)*PP;
        A := [A[1:stairStep, 1:stairStep],A[1:stairStep, stairStep + 1:nx]*Q2;
          transpose(Q2)*A[stairStep + 1:nx, 1:stairStep],transpose(Q2)*A[
          stairStep + 1:nx, stairStep + 1:nx]*Q2];
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
      stairStep := if Modelica.Math.Matrices.isEqual(ss.B, zeros(size(ss.B, 1), size(ss.B, 2))) then 0 else 
              1;
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
    else // no inputs, nu==0
    isControllable := false;
    ssm1 :=  Modelica_LinearSystems2.Internal.StateSpaceR(
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
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    output Boolean isControllable;
    output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
      redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "upper block Hessenberg form state space system";
    output Real P[:,:];

    protected
    Real A[size(ss.A, 1),size(ss.A, 2)];
    Real B[size(ss.B, 1),size(ss.B, 2)];
    Real C[size(ss.C, 1),size(ss.C, 2)];
    Real U[:,:];
    Real VT[:,:];

    Real Q[:,:];
    Real Pi[:,:];
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

      B := [diagonal(sigma[1:rankS]),zeros(rankS, nu - rankS); zeros(nx - rankS, nu)];
      P := transpose(U);
      Q := transpose(VT);
      A := transpose(U)*ss.A*U;
  //    C := ss.C*U;
      (sigma,U,VT) := Modelica.Math.Matrices.singularValues(A[rankS + 1:nx, 1:rankS]);

      stairStep := rankS;
      rankS := 0;
      if size(sigma, 1) > 1 then
        for i in 1:size(sigma, 1) loop
          if sigma[i] > maxA*sigma[1]*Modelica.Constants.eps then
            rankS := rankS + 1;
          end if;
        end for;
      else
        rankS := if size(sigma, 1) > 0 then if sigma[1] > maxA*Modelica.Constants.eps then 1 else 0 else 0;
      end if;

      stairStep := stairStep + rankS;
      Pi := [VT,zeros(size(VT, 1), size(U, 1)); zeros(size(U, 2), size(VT, 2)),transpose(U)];
      B := Pi*B;
      P := Pi*P;

   // P could be ambigious according to the sign
      if transpose(P)*B[:,1]*ss.B[:,1]<0 then
        Pi:=-Pi;
        P:=-P;
      end if;

      A := Pi*A*transpose(Pi);

     // should be made better because of many zeros in B

      while stairStep < nx and rankS > 0 and not 
          Modelica.Math.Matrices.isEqual(
              A[stairStep + 1:nx, stairStep - rankS + 1:stairStep],
              zeros(stairStep - rankS, rankS),
              eps) loop

        (sigma,U,VT) := Modelica.Math.Matrices.singularValues(A[stairStep + 1:nx, stairStep
           - rankS + 1:stairStep]);

        Pi := [identity(stairStep - rankS),zeros(stairStep - rankS, nx -
          stairStep + rankS); zeros(nx - stairStep + rankS, stairStep - rankS),
          [VT,zeros(rankS, nx - stairStep); zeros(nx - stairStep, rankS),
          transpose(U)]];
        P := Pi*P;
        A := Pi*A*transpose(Pi);

  //new implenmentation advisable because of many zeros in Pi
  //      C := C*transpose(Pi);
        rankS := 0;
        if size(sigma, 1) > 1 then
          for i in 1:size(sigma, 1) loop
            if sigma[i] > maxA*sigma[1]*Modelica.Constants.eps then
              rankS := rankS + 1;
            end if;
          end for;
        else
          rankS := if size(sigma, 1) > 0 then if sigma[1] > maxA*Modelica.Constants.eps then 
                  1 else 0 else 0;
        end if;
        stairStep := stairStep + rankS;
      end while;

      B := P*ss.B;
      C := ss.C*transpose(P);

    else
      stairStep := if Modelica.Math.Matrices.isEqual(ss.B, zeros(size(ss.B, 1),
        size(ss.B, 2))) then 0 else 1;
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

    else // no inputs, nu==0
    isControllable := false;
    ssm1 :=  Modelica_LinearSystems2.Internal.StateSpaceR(
          A=ss.A,
          B=ss.B,
          C=ss.C,
          D=ss.D,
          r=0);
    P := identity(nu);
     end if;

  end staircaseSVD;

  encapsulated function partialGain
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.StateSpace.Internal;

    input Real H[:,size(H, 1)];
    input Real b[size(H, 1)];
    output Real result;
    protected
    Real Hh[:,:]=H;
    Real bh[:]=b;
    Integer q=size(H, 1);
  algorithm

    (Hh,bh) := Internal.trianUpperHess(Hh, bh);
    result := bh[q]/Hh[q, q];

  end partialGain;

  encapsulated function assignOneOrTwoPoles_alpha
      "Algorithm to assign p (p = 1 or 2) eigenvalues"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Math.Vectors;

    input Real F[:,size(F, 1)] "system matrix of order p=1 or p=2";
    input Real G[size(F, 1),:] "control input matrix p rows";
    input Complex gamma[size(F, 1)];
    input Real tolerance=Modelica.Constants.eps;
    output Real K[:,size(F, 1)] "feedback matrix p columns";

    protected
    Real Gamma[:,:];
    Integer rankGs;
    Real Fs[size(F, 1),size(F, 2)];
    Real Gs[size(G, 1),size(G, 2)];
    Real Gst[:,:]=transpose(G);
    Real Ks[:,size(F, 1)];
    Real c;
    Real s;
    Real r;

    Real V1[size(G, 2),size(G, 2)];
    Real V2[size(G, 2),size(G, 2)];
    Real V[size(G, 2),size(G, 2)];
    Real U[size(F, 1),size(F, 2)];

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
      assert(gamma[1].im == 0, "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F has size [" + String(size(F, 1)) + "," + String(size(F, 1)) +
        "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
         + String(gamma[1].im));
    elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
      assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im,
        "\nThe assigned pole pair given in function StateSpace.Internal.assignOneOrTwoPoles() must be conjungated complex. However, the poles are\npole1 = "
         + String(gamma[1]) + "\npole2 = " + String(gamma[2]) +
        ". \nTry\npole1 = " + String(gamma[1]) + "\npole2 = " + String(
        Complex.conj(gamma[1])) + "\ninstead");
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
          U := [c,s; -s,c];
        end if;
        Gs := U*G;

        rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
      else
       // size(G, 2)>1

        if size(G, 1) == 1 then // U=I, compute V by just one Householder transformation
          U := [1];
          u1 := cat(1, Vectors.householderVector(Gst[:, 1], cat(1, {1}, zeros(size(G, 2) - 1))));
                                   // Householder vector
          Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(
            Gst, u1);

          V := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1*u1);
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
          if Modelica.Math.Vectors.isEqual(Gst[:, 2], zeros(size(G, 2)), tolerance) or 
            Modelica.Math.Matrices.isEqual(Gst[2:size(Gst, 1), :], zeros(size(Gst, 1) - 1, size(Gst, 2)), tolerance) then
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
          U := [c,s; -s,c];
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
          Ks := [(Fs[1, 1] - gamma[1].re)/Gs[1, 1] - Gs[1, 2]*(Fs[2, 1] - 1)/Gs[1, 1]/Gs[2,
            2],(Fs[1, 2] + (gamma[1].im)^2)/Gs[1, 1] - Gs[1, 2]*(Fs[2, 2] -
            gamma[2].re)/Gs[1, 1]/Gs[2, 2]; (Fs[2, 1] - 1)/Gs[2, 2],(Fs[2, 2]
             - gamma[2].re)/Gs[2, 2]];
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

  encapsulated function assignOneOrTwoPoles
      "Algorithm to assign p (p = 1 or 2) eigenvalues"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Math.Vectors;

    input Real F[:,size(F, 1)] "system matrix of order p=1 or p=2";
    input Real G[size(F, 1),:] "control input matrix p rows";
    input Complex gamma[size(F, 1)];
    input Real tolerance=Modelica.Constants.eps;
    output Real K[:,size(F, 1)] "feedback matrix p columns";

    protected
    Real Gamma[:,:];
    Integer rankGs;
    Real Fs[size(F, 1),size(F, 2)];
    Real Gs[size(G, 1),size(G, 2)];
    Real Gst[:,:]=transpose(G);
    Real Ks[:,size(F, 1)];
    Real c;
    Real s;
    Real r;
    Integer p=size(F,1);
    Real sigmaG[:];

    Real V[size(G, 2),size(G, 2)];
    Real U[size(F, 1),size(F, 2)];

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
      assert(gamma[1].im == 0, "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F has size [" + String(size(F, 1)) + "," + String(size(F, 1)) +
        "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
         + String(gamma[1].im));
    elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
      assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im,
        "\nThe assigned pole pair given in function StateSpace.Internal.assignOneOrTwoPoles() must be conjungated complex. However, the poles are\npole1 = "
         + String(gamma[1]) + "\npole2 = " + String(gamma[2]) +
        ". \nTry\npole1 = " + String(gamma[1]) + "\npole2 = " + String(
        Complex.conj(gamma[1])) + "\ninstead");
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
          U := [c,s; -s,c];
        end if;
        Gs := U*G;

        rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
      else
       // size(G, 2)>1

        if size(G, 1) == 1 then // U=I, compute V by just one Householder transformation
          U := [1];
          u1 := cat(1, Vectors.householderVector(Gst[:, 1],
                       cat(1, {1}, zeros(size(G, 2) - 1))));// Householder vector
          Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(Gst, u1);

          V := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1*u1);
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
        tolerance), "A subsystem in StateSpace.Internal.assignOneOrTwoPoles() is not controllable");

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

          Ks := [(Fs[1, 1] - gamma[1].re)/Gs[1, 1],  (Fs[1, 2] + (gamma[1].im)^2)/Gs[1, 1];
          (Fs[2, 1] - 1)/Gs[2, 2],(Fs[2, 2] - gamma[2].re)/Gs[2, 2]];
        else

          Ks[1, 1] := (gamma[1].re + gamma[2].re - Fs[1, 1] - Fs[2, 2])/Gs[1, 1];
          Ks[1, 2] := Ks[1, 1]*Fs[2, 2]/Fs[2, 1] + (Fs[1, 1]*Fs[2, 2] - Fs[1, 2]*
            Fs[2, 1] - (gamma[1].re*gamma[2].re - gamma[1].im*gamma[2].im))/Fs[2,1]/Gs[1, 1];
          Ks := -Ks;
        end if;
      end if;

      K := V[:, 1:size(Ks, 1)]*Ks*U;

    else
      if p == 1 then
        Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, "   + String(F[1, 1]) + " remains and " + String(gamma[1].re) + " cannot be realized");
      else
        system_ev := Complex.eigenValues(F);
        Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, "   + String(system_ev[1].re) + (if abs(system_ev[1].im) > 0 then " + " else 
                " - ") + String(system_ev[1].im) + "j and " + String(system_ev[2].re)
           + (if abs(system_ev[2].im) > 0 then " + " else " - ") + String(
          system_ev[2].im) + "j remain and " + String(gamma[1].re) + (if abs(
          gamma[1].im) > 0 then (if gamma[1].im > 0 then " + " else " - " +
          String(gamma[1].im) + "j") else "" + " and ") + String(gamma[2].re) + (
          if abs(gamma[2].im) > 0 then (if gamma[2].im > 0 then " + " else " - " +
          String(gamma[2].im) + "j") else "") + " cannot be realized");
      end if;
      K := zeros(size(G, 2), size(F, 1));
    end if;

  end assignOneOrTwoPoles;

  encapsulated function readSystemDimension
      "Read the order nx of state matrix and the numbers nu and ny of inputs and outputs"
      import Modelica_LinearSystems2;
    input String fileName="stateSpace.mat" 
                                annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    input String matrixName="ABCD"
        "Name of the generalized state space system matrix";
    output Integer xuy[3];

    protected
    Real sizeA[1,1]=readMatrix(
            fileName,
            "nx",
            1,
            1);

    Integer ABCDsizes[2]=readMatrixSize(fileName, matrixName);

  algorithm
    xuy[1] := integer(sizeA[1, 1]);
    xuy[2] := ABCDsizes[2] - xuy[1];
    xuy[3] := ABCDsizes[1] - xuy[1];

  end readSystemDimension;

  encapsulated function readLength_nu
      "Read the number of inputs nu of a state space system from a file"

    input String fileName="ss_siso.mat"
        "Name of the state space system data file"     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    input String matrixName="ABCD" "Name of the state space system matrix"    annotation(Dialog);

    output Integer nu;
    protected
    Real nxMat[1,1]=readMatrix(
            fileName,
            "nx",
            1,
            1);
    Integer ABCDsizes[2]=readMatrixSize(fileName, matrixName);
    Integer nx=integer(nxMat[1, 1]);

  algorithm
    nu := ABCDsizes[2] - nx;
  end readLength_nu;

  encapsulated function readLength_ny
      "Read the number of outputs ny of a state space system from a file"

    input String fileName="ss_siso.mat"
        "Name of the state space system data file"     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    input String matrixName="ABCD" "Name of the state space system matrix"    annotation(Dialog);

    output Integer ny;
    protected
    Real nxMat[1,1]=readMatrix(
            fileName,
            "nx",
            1,
            1);
    Integer ABCDsizes[2]=readMatrixSize(fileName, matrixName);
    Integer nx=integer(nxMat[1, 1]);

  algorithm
    ny := ABCDsizes[1] - nx;
  end readLength_ny;

  encapsulated function dgreeOfRedSys
      "Calculate the controllable and observable part of a state space system"

      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2;

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

  encapsulated function householder
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

    input StateSpace ss;
    input Real v[size(ss.A, 1)];

    output StateSpace ssh;

    protected
    Real Ah[size(ss.A, 1),size(ss.A, 2)];
    Real Bh[size(ss.B, 1),size(ss.B, 2)];
    Real Ch[size(ss.C, 1),size(ss.C, 2)];
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

  encapsulated function complexZeros
      "Calculate the zeros of the related transfer function"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss "StateSpace object";
    output Complex zeros2[:]=fill(Complex(0, 0),
        StateSpace.Internal.numberOfZeros(ss));

    protected
    Integer nx=size(ss.A, 2);

    Complex zeros[:];
    Complex poles[:];
    Real eval[nx,2];
    Real evec[nx,nx];

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

    zeros := StateSpace.Analysis.invariantZeros( ss);

    poles := Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A);

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

  encapsulated function numberOfPoles
      "Calculate the number of poles of the related transfer function"

      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss "StateSpace object";

    output Integer numberOfPoles=StateSpace.Internal.numberOfPolesAndZeros(ss);

  algorithm
  end numberOfPoles;

  encapsulated function numberOfPolesAndZeros
      "Calculate the number poles and of zeros of the related transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

    input StateSpace ss "StateSpace object";
    output Integer numberOfPoles;
    output Integer numberOfZeros;

    protected
    Integer nx=size(ss.A, 2);
    Complex zeros[:];
    Complex zeros2[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";
    Complex poles[:];
    Complex poles2[:] "eigenvalues of ss";
    Real eval[nx,2];
    Real evec[nx,nx];

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

    poles := Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A);
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

  encapsulated function complexPoles
      "Generate a zeros-and-poles representation from state space representation"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss "StateSpace object";

    output Complex poles2[:]=fill(Complex(0, 0),
        StateSpace.Internal.numberOfPoles(ss));

    protected
    Integer nx=size(ss.A, 2);

    Complex zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";
    Complex zeros2[:];
    Complex poles[:] "eigenvalues of ss";
    Real eval[nx,2];
    Real evec[nx,nx];

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
    zeros := StateSpace.Analysis.invariantZeros( ss);
    zeros2 := zeros;

    poles := Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A);

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

  encapsulated function trianUpperHess
      "Triangulize an upper Hessenberg matrix by repeatedly applicated householder reflexion"
      import Modelica;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

    input Real H[:,:];
    input Real b[size(H, 1)];

    output Real Ht[size(H, 1),size(H, 2)];
    output Real bt[size(b, 1)];

    protected
    Integer q=size(H, 1);
    Real u[:] "householder vector";
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
This function computes an triangular matrix from an upper Hessenberg matrix by stepwise annihilation of the subdiagonal elements.
<p>
<blockquote><pre>
<b>A</b> -> <b>QA</b> = <b>T</b>
</pre></blockquote>
</p>
It is assumend that the original matrix has upper hessenberg form.
Additionally the vector b is transformed in the same way
<blockquote><pre>
<b>b</b> -> <b>Qb</b> = <b>q</b>
</pre></blockquote>
</p>
The function is primarily used to calculate the transfer function gain from a SISO state space system in observer Hessenberg form
<blockquote><pre>
     ( *   *   ...   ...    * )          ( * )
     ( *   *   ...   ...    * )          ( . )
 <b>A</b> = ( 0   *   ...   ...    * ),    <b>b</b> =  ( . ),   <b>c</b> = ( 0, ..., 0, * )
     ( .   .    .     .     . )          ( * )
     ( 0  ...   0     *     * )          ( * )
 
</pre>
</blockquote>
If <b>A</b> is upper Hessenberg and <b>T</b> = <b>Q</b>*<b>A</b> is triangular then obviously <b>H</b>(s) = <b>Q</b>*(s*<b>I</b> -<b>A</b>) = s*<b>I</b> - <b>T</b>.
<p>
Further on, if <b>T</b> is triangular then also <b>H</b> = s<b>I</b> - <b>T</b> is and the element l_nn of <b>L</b> = inv(<b>H</b>) is given by 1/h_nn.
The frequency response G(s0)for a given s0 that is neither zero nor pole of the system can be calculated by
<blockquote><pre>
                    -1               -1     -1         -1           -1  -1        -1
G(s0)  = <b>c</b>*(s0*<b>I</b> -<b>A</b>)  *<b>b</b> = <b>c</b>*(s0*<b>I</b> -<b>A</b>)  *<b>Q</b>*<b>Q</b>  *<b>b</b> = <b>c</b>*(<b>Q</b>  *(s0*<b>I</b> -<b>A</b>))  *<b>Q</b>  *<b>b</b> = <b>c</b>*<b>H</b>(s0)*<b>q</b>
</pre></blockquote>
and because only the n'th element of <b>c</b> is different to zero the gain k is given by
<blockquote><pre>
    q_nn*c_nn     product(s0 - poles_i)
k = ---------- * ----------------------
       h_nn       product(s0 - zeros_i)
</pre></blockquote>
 
 
</p>
 
</html>"));
  end trianUpperHess;

  encapsulated function reducedCtrSystem2
      "calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;

    input StateSpace ss;
    input Real eps=0;

    protected
    StateSpace sst(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
      redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "tranformed state space system";

    Integer nx=size(ss.A, 1);
    Real Ah1[size(ss.A, 1),size(ss.A, 2)];
    Real bh1[size(ss.A, 1)];
    Real ch1[size(ss.A, 1)];

    Real u[:] "householder vector";
    Real cpoles[:,2]=Modelica_LinearSystems2.StateSpace.Internal.controllablePoles(ss);
    Integer rankQc=size(cpoles,1);

    Integer rankQc2;
    Real Qc2[nx,nx];
    Real sigma[:];
    Real eps2;

    Real Ah2[rankQc,rankQc];
    Real bh2[rankQc];
    Real ch2[rankQc];

    Integer ll;
    Integer r;

    Boolean h;

    public
    output StateSpace ssm1(
      redeclare Real A[rankQc,rankQc],
      redeclare Real B[rankQc,1],
      redeclare Real C[1,rankQc],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
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
        D=[0]);

    output Real cPoles[:,2] "controllable poles";
    output Real ncPoles[:,2] "uncontrollable poles";
    output Real poles[size(ss.A, 1),2] "controllable and uncontrollable poles";
    protected
    Modelica_LinearSystems2.Internal.StateSpaceR ssch(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
      redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "upper block controller Hessenberg form state space system";
    Boolean isControllable;

  algorithm
    if size(ss.B, 2) == 0 then
      poles := Modelica.Math.Matrices.eigenValues(ss.A);
      ncPoles := poles;
      cPoles := fill(0, 0, 2);
    else
  // build upper Hessenberg staircase to decomposite controllable/uncontrollable subspaces
  // The controllable part of A is in A[1:ssch.r, 1:ssch.r]
      (isControllable,ssch) := StateSpace.Internal.staircaseSVD(ss);
      if isControllable then
        poles := Modelica.Math.Matrices.eigenValues(ss.A);
        cPoles := poles;
        ncPoles := fill(0, 0, 2);
      else
        cPoles := Modelica.Math.Matrices.eigenValues(ssch.A[1:ssch.r, 1:ssch.r]);
        ncPoles := Modelica.Math.Matrices.eigenValues(ssch.A[ssch.r + 1:size(ss.A,
          1), ssch.r + 1:size(ss.A, 1)]);
        poles := [cPoles; ncPoles];
      end if;
    end if;

    annotation (Documentation(info="<html>
The function uses the SVD based staircase algorithm to transform the state space representation into a similar state space
to separate the uncontrollable poles from the controllable poles.
</html>"));
  end controllablePoles;

  encapsulated function polesAndZeros
      "Generate poles and zeros from state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Internal.PolesAndZeros;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Internal;

    input StateSpace ss "StateSpace object";
    input StateSpace ssm=Modelica_LinearSystems2.StateSpace.Transformation.toIrreducibleForm(
                                                                         ss);
    output Internal.PolesAndZeros pz(
      redeclare Real p_real[size(ssm.A, 1)],
      redeclare Real p_im[size(ssm.A, 1)],
      redeclare Real z_real[size(StateSpace.Analysis.invariantZeros(ssm), 1)],
      redeclare Real z_im[size(StateSpace.Analysis.invariantZeros(ssm), 1)]);
    protected
    Complex poles[:]=Complex.eigenValues(ssm.A);
    Complex zeros[:]=StateSpace.Analysis.invariantZeros( ssm);

  algorithm
    pz.p_real := poles[:].re;
    pz.p_im := poles[:].im;
    pz.z_real := zeros[:].re;
    pz.z_im := zeros[:].im;
    pz.norz_p := Internal.numberOfRealZeros(poles);
    pz.norz_z := Internal.numberOfRealZeros(zeros);

  end polesAndZeros;

  encapsulated function scaleFactor1
      "Return scale factor for first order block"
      import Modelica;
    input Real n "(s+n)/(s+d)";
    input Real d "(s+n)/(s+d)";
    input Real small=100*Modelica.Constants.eps;
    output Real k "= d/n, if d,n are not zero, otherwise special cases";
  algorithm
  //  k := (if abs(d) > small then abs(d) else 1)/(if abs(n) > small then abs(n) else 1);
    k := if abs(d) > small  and abs(n) > small then abs(d)/abs(n) else 1;

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

  encapsulated function invariantZerosHessenberg
      "Fast version to calculate the system zeros of a SISO system with D=[0] and A has upper Hessenberg form, delivered by StateSpace.reduceSystem"
      import Modelica;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

    input StateSpace ss "Linear system in state space form";
    output Complex Zeros[:]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
    Integer nx=size(ss.A, 1) "Number of states";
    Integer nu=size(ss.B, 2) "Number of inputs";
    Integer ny=size(ss.C, 1) "Number of outputs";
    Real Ah[nx,nx]=ss.A;
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
        Zeros := Complex.Internal.eigenValues_dhseqr(Ah[1:k - 1, 1:k - 1]);

        for i in 1:k - 1 loop
          if Complex.'abs'(Zeros[i]) < Modelica.Math.Matrices.norm(Ah[1:k - 1, 1:
              k - 1], p=1)*1e-12 then
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

  encapsulated function cntrHessenberg
      "calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),1],
      redeclare Real C[1,size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "controllable state space system";

    protected
    Integer nx=size(ss.A, 1);
    Real Ah1[nx,nx];
    Real bh1[nx];
    Real ch1[nx];
    Real u[:] "householder vector";
    Real Q[nx,nx];
    Real V[size(ss.A, 1),size(ss.A, 2)];
    Real tau[nx - 1];
    Real Qc[:,:];
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
              fill(0, nx - 1)));  //householder vector to compute a housholder reflector S = I - 2*u*u'/u'*u
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

  encapsulated function transposeStateSpace
      "Return the transposed state space system"

      import Modelica_LinearSystems2.StateSpace;

    input StateSpace ss;

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

    input Real A[:,:];
    input Real B[:,:];
    input Real C[:,:];
    input Real D[:,:];

    output Real Ar[:,:];
    output Real Br[:,:];
    output Real Cr[:,:];
    output Real Dr[:,:];
    output Integer n;
    output Integer m;
    output Integer p;

    protected
    Real A2[:,:];
    Real B2[:,:];
    Real C2[:,:];
    Real CC[:,:];
    Real Co[:,:];
    Real Cu[:,:];
    Real D2[:,:];
    Real DD[:,:];
    Real Mr[:,:];
    Real Vf[:,:];
    Real V[:,:];
    Real V2[:,:];
    Real R[:,:];
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
    Real normA=Modelica.Math.Matrices.norm(A=A, p=1);
    Real eps=normA*1e-12;

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
        (V,R,tau,V2) := Matrices.QR( D2);

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
        stop1 := size(CC,1) == rankR;

        if not stop1 then
          Cu := CC[sigma + 1:end, :];
          Co := CC[1:sigma, :];

          (V,R,tau,V2) := Matrices.QR( Matrices.fliplr(transpose(Cu)));
           Vf:=Matrices.fliplr(V2);

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
          stop2 := size(Cu,2) == rankR;

          if not stop1 and not stop2 then
            nue := size(Cu, 2) - rankR;
            mue := rho + sigma;
            delta := delta + rho;

            if sigma == 0 then
              Mr := [transpose(Vf)*A2*Vf,transpose(Vf)*B2];
            else
              Mr := [transpose(Vf)*A2*Vf,transpose(Vf)*B2; Co*Vf,DD];
            end if;

            A2 := Mr[1:nue, 1:nue];
            B2 := Mr[1:nue, nue + rho + 1:nue + rho + nu];
            C2 := Mr[nue + 1:nue + mue,1:nue];
            D2 := Mr[nue + 1:nue + mue,nue + rho + 1:nue + rho + nu];

            j := j + 1;
          end if;
         //not stop1 or not stop2

        end if;
         //if not stop1

        stop := stop1 or stop2 or j>3*nx;

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
        Ar := fill(0,0,0);
        Br := fill(0,0,0);
        Cr := fill(0,0,0);
        Dr := fill(0,0,0);
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

    annotation (Documentation(info="<html></html>"));
  end reduceRosenbrock;

  encapsulated function reducedCtrSystem
      "calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),1],
      redeclare Real C[1,size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "controllable state space system";

    protected
    Integer nx=size(ss.A, 1);
    Real Ah1[nx,size(ss.A, 2)];
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

        while r <= nx - 1 and maxa > nx*normA*1e-10 loop
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
      "calculate the controllable part of a SISO system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math.Complex;

    input StateSpace ss;

    output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
      redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
      redeclare Real B[size(ss.B, 1),1],
      redeclare Real C[1,size(ss.C, 2)],
      redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
        "controllable state space system";

    protected
    Integer nx=size(ss.A, 1);
    Real Ah1[nx,size(ss.A, 2)];
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

end Internal;

    annotation (
    defaultComponentName="stateSpace",
    Documentation(info="<html>
<p>
This record defines a linear time invariant differential
equation system in state space form:
</p>
<pre>    <b>der</b>(x) = A * x + B * u
        y  = C * x + D * u
</pre>
<p>
with
</p>
<ul>
<li> u - the input vector</li>
<li> y - the output vector</li>
<li> x - the state vector</li>
<li> A,B,C,D - matrices of appropriate dimensions</li>
</ul>
</html>"),
    DymolaStoredErrors);
end StateSpace;
