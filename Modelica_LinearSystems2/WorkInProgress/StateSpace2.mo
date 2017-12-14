within Modelica_LinearSystems2.WorkInProgress;
record StateSpace2
  "Continuous state space description of a linear, time invariant differential equation system (data + operations)"

  extends Modelica.Icons.Record;

  Real A[:,size(A, 1)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real B[size(A, 1),:]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real C[:,size(A, 1)]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));

//   String uNames[size(B, 2)]=fill("", size(B, 2))  annotation(Dialog(group="Signal names"));
//   String yNames[size(C, 1)]=fill("", size(C, 1)) annotation(Dialog(group="Signal names"));
//   String xNames[size(A, 1)]=fill("", size(A, 1))  annotation(Dialog(group="Signal names"));

encapsulated operator 'constructor'
    "Default constructors for a StateSpace record"
    import Modelica_LinearSystems2;

  function fromABCDMatrices "Default constructor for a StateSpace record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.StateSpace2;

    input Real A[:,size(A, 1)];
    input Real B[size(A, 1),:];
    input Real C[:,size(A, 1)];
    input Real D[size(C, 1),size(B, 2)];

      //     input String uNames[size(B, 2)]=fill("", size(B, 2));
      //     input String yNames[size(C, 1)]=fill("", size(C, 1));
      //     input String xNames[size(A, 2)]=fill("", size(A, 2));

     output StateSpace2 result(
      redeclare Real A[size(A, 1),size(A, 2)],
      redeclare Real B[size(B, 1),size(B, 2)],
      redeclare Real C[size(C, 1),size(C, 2)],
      redeclare Real D[size(D, 1),size(D, 2)]);

      //       redeclare String uNames[size(B, 2)],
      //       redeclare String yNames[size(C, 1)],
      //       redeclare String xNames[size(A, 2)]);

  algorithm
    result.A := A;
    result.B := B;
    result.C := C;
    result.D := D;
      //      result.uNames := uNames;
      //      result.yNames := yNames;
      //      result.xNames := xNames;

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  ss </td><td align=center>=</td>  <td> 'constructor'.<b>fromABCDMatrices</b>(A, B, C, D)  </td> </tr>
</table>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss with<br>
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

  function fromReal "Generate a StateSpace data record from a Real value"

      import Modelica;
      import Modelica_LinearSystems2.WorkInProgress.StateSpace2;

    input Real r "Value of Real variable";
    output StateSpace2 ss(
      redeclare Real A[0,0],
      redeclare Real B[0,1],
      redeclare Real C[1,0],
      redeclare Real D[1,1]) "= r";

  algorithm
    ss.D[1, 1] := r;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  ss </td><td align=center>=</td>  <td> 'constructor'.<b>fromReal</b>(r)  </td> </tr>
</table>

<h4>Description</h4>
<p>
This function constructs a StateSpace record ss from a Real value, i.e. a state space system without a state and an output without dynamics:
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

// algorithm
//            // this is the constructor algorithm
//   result.A := A;
//   result.B := B;
//   result.C := C;
//   result.D := D;
//   //end constructor;

  function fromTransferFunction =
      Modelica_LinearSystems2.TransferFunction.Conversion.toStateSpace;
  function fromZerosAndPoles =
      Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace;

    annotation (Documentation(info="<html>
This package contains the default constructors for StateSpace record.
</html>"));
end 'constructor';

encapsulated operator '-'
    "Contains operators for subtraction of state space systems"

  function subtract
      "Subtraction of two state space systems connected in parallel (= inputs are the same, outputs of the two systems are subtracted)"

      import Modelica;
      import Modelica_LinearSystems2.WorkInProgress.StateSpace2;

    input StateSpace2 ss1 "State space system 1";
    input StateSpace2 ss2 "State Space system 2 is subtracted from system 1";
    output StateSpace2 result(
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
<h4>Syntax</h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> Modelica_LinearSystems2.StateSpace.'-'.<b>subtract</b>(ss1, ss2)  </td> </tr>
</table>

<h4>Description</h4>
<p>
This operator function computes the subtraction of two state space systems connected in parallel, i.e. the inputs are the same and the outputs of the two systems are subtracted. Therefore, The systems must have the same number of inputs and outputs but not the same number of states. The resulting system has an order of system_order1 + system_order2.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
ss3 := ss1 - ss2;
</pre></blockquote>

<h4>Example</h4>
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

</html>"));
  end subtract;

  function negate
      "Unary minus (state space system where the output is multiplied by a gain of -1)"
      import Modelica;
      import StateSpace = Modelica_LinearSystems2.WorkInProgress.StateSpace2;

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
<h4>Description</h4>
<p>
This package contains the <a href=\"modelica://Modelica_LinearSystems2.StateSpace.'-'.subtract#info\">'subtract'</a> and the <a href=\"modelica://Modelica_LinearSystems2.StateSpace.'-'.subtract#info\">'-'</a> operator for StateSpace records.

</html>"));
end '-';

encapsulated operator function '+'
    "Parallel connection of two state space systems (= inputs are the same, outputs of the two systems are added)"
    import StateSpace = Modelica_LinearSystems2.WorkInProgress.StateSpace2;

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
  //end '+';
//  end add;

end '+';

encapsulated operator function '*'
    "Series connection of two state space systems"
    import StateSpace = Modelica_LinearSystems2.WorkInProgress.StateSpace2;

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
//  end multiply;

end '*';

encapsulated operator function '=='
    "Check whether two linear systems have identical matrices"
    import Modelica.Math.Matrices.isEqual;
    import StateSpace = Modelica_LinearSystems2.WorkInProgress.StateSpace2;

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

encapsulated package Import
    "Utilitiy functions to import StaeSpace representations"

  encapsulated function fromFile "Read a StateSpace data record from mat-file"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.WorkInProgress.StateSpace2;
    import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

    input String fileName="dslin.mat"
        "Name of the state space system data file"     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));
    input String matrixName="ABCD" "Name of the state space system matrix"    annotation(Dialog);
    protected
    input Integer xuy[3] = ReadSystemDimension(fileName, matrixName);
    input Integer nx=xuy[1];
    input Integer nu=xuy[2];
    input Integer ny=xuy[3];

    public
    output StateSpace2 result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "= model linearized at initial point";

    protected
    Real ABCD[nx + ny,nx + nu] = Modelica.Utilities.Streams.readRealMatrix(
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
<h4>Syntax</h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> StateSpace.Import.<b>fromFile</b>(fileName, matrixName)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Reads and loads a state space system from a mat-file <tt>fileName</tt>. The file must contain the matrix [A, B; C, D] named matrixName and the integer nx representing the order of the system, i.e. the number of rows of the square matrix A.

<h4>Example</h4>
<blockquote><pre>
<b>algorithm</b>
  ss:=Modelica_LinearSystems2.StateSpace.Import.fromFile(&quot;stateSpace.mat&quot;, &quot;ABCD&quot;);
//  ss=StateSpace(
      A=[-1, 0, 0; 0, -2, 0; 0, 0, -3],
      B=[1; 1; 0],
      C=[1, 1, 1],
      D=[0])


</pre></blockquote>


</html>"));
  end fromFile;

    annotation (Documentation(info="<html>
</html>"));

end Import;

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
</html>"));
end StateSpace2;
