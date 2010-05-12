within Modelica_LinearSystems2;
record DiscreteStateSpace
  "Discrete state space description of a linear, time invariant difference equation system (data + operations)"
  extends Modelica.Icons.Record;

  Real A[:,size(A, 1)]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real B[size(A, 1),:]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real C[:,size(A, 1)]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));

  Modelica.SIunits.Time Ts=1 "Sample time" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));
  Real B2[size(B, 1),size(B, 2)]=fill(0,size(B,1),size(B,2))
    "Reconstruct continuous state" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));
  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method" 
        annotation(Dialog(group="Data used to construct discrete from continuous system"));

//      String yNames[size(C, 1)]=fill("", size(C, 1)) "Names of the output signals" annotation(Dialog(group="Signal names"));
//      String xNames[size(A, 1)]=fill("", size(A, 1)) "Names of the states"  annotation(Dialog(group="Signal names"));
//      String uNames[size(B, 2)]=fill("", size(B, 2)) "Names of the input signals" annotation(Dialog(group="Signal names"));

  encapsulated operator 'constructor'
    "Default constructors for a DiscreteStateSpace record"
  import Modelica_LinearSystems2;
  function fromDiscreteTransferFunction = 
      Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteStateSpace
                                                                                             annotation (Documentation(info="<html> </html>"));
  function fromDiscreteZerosAndPoles = 
        Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteStateSpace
                                                                                                          annotation (Documentation(info="<html> </html>"));

     function fromReal "Generate a StateSpace data record from a Real value"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

        input Real r "Value of Real variable";
        input Modelica.SIunits.Time Ts=1 "Sample time";
        input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
        output DiscreteStateSpace dss(
      redeclare Real A[0,0],
      redeclare Real B[0,1],
      redeclare Real B2[0,1],
      redeclare Real C[1,0],
      redeclare Real D[1,1]) "= r";

     algorithm
        dss.D[1, 1] := r;
        dss.Ts := Ts;
        dss.method := method;
        annotation (overloadsConstructor=true, Documentation(info=
                                                              "<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dss </td><td align=center>=</td>  <td> 'constructor'.<b>fromReal</b>(r)  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function constructs a DiscreteStateSpace record dss from a Real value, i.e. a discrete state space system without a state and an output without dynamics:
<blockquote><pre>
y = r*u
</pre></blockquote>
Therefore, the matrices are defined by
<blockquote><pre>
  dss.A = fill(0,0,0);
  dss.B = fill(0,0,1);
  dss.C = fill(0,1,0);
  ss.D = [r];
  dss.B2 = fill(0,0,1);
</pre></blockquote>
 
</p>
 
 
</html>"));
     end fromReal;

    function fromMatrices "Default constructor for a DiscreteStateSpace record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input Real A[:,size(A, 1)] 
                              annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real B[size(A, 1),:] 
                              annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real C[:,size(A, 1)] 
                              annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real D[size(C, 1),size(B, 2)] 
                                      annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));

      input Modelica.SIunits.Time Ts=1 "Sample time" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));
      input Real B2[:,:]=zeros(size(B, 1), size(B, 2));
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization methodDiscretization method" 
        annotation(Dialog(group="Data used to construct discrete from continuous system"));
      output DiscreteStateSpace result(
        redeclare Real A[size(A, 1),size(A, 2)],
        redeclare Real B[size(B, 1),size(B, 2)],
        redeclare Real C[size(C, 1),size(C, 2)],
        redeclare Real D[size(D, 1),size(D, 2)],
        redeclare Real B2[size(B2, 1),size(B2, 2)]);
    algorithm
      result.A := A;
      result.B := B;
      result.C := C;
      result.D := D;

      result.B2 := B2;
      result.Ts := Ts;
      result.method := method;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right> dss </td><td align=center>=</td>  <td> 'constructor'.<b>fromMatrices</b>(A, B, C, D)  </td> </tr>
<tr> <td align=right> dss </td><td align=center>=</td>  <td> 'constructor'.<b>fromMatrices</b>(A, B, C, D, Ts, B2, method)  </td> </tr>

</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function constructs a DiscreteStateSpace record dss with<br>
<blockquote><pre>
  dss.A = A;
  dss.B = B;
  dss.C = C;
  dss.D = D;
  dss.B2 = B2;
  dss.Ts = Ts;
  dss.method = method;
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
  dss := 'constructor'.fromMatrices(A, B, C, D);
  // dss.A = [1]
  // dss.B = [1]
  // dss.C = [1]
  // dss.D = [0]
  // dss.B2 = [0]
  // dss.Ts = 1
  // dss.method = Modelica_LinearSystems2.Types.Method.Trapezoidal 
  
  

</pre></blockquote>
</html>"));
    end fromMatrices;

    function fromStateSpace
      "Transform a continuous into a discrete linear state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Types.Method;
      import Modelica_LinearSystems2.Math.Matrices.LU_solve2;

      input Modelica_LinearSystems2.StateSpace sc
        "Continuous linear state space system";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      output Modelica_LinearSystems2.DiscreteStateSpace sd(
        redeclare Real A[size(sc.A, 1),size(sc.A, 2)],
        redeclare Real B[size(sc.B, 1),size(sc.B, 2)],
        redeclare Real C[size(sc.C, 1),size(sc.C, 2)],
        redeclare Real D[size(sc.D, 1),size(sc.D, 2)],
        redeclare Real B2[size(sc.B, 1),size(sc.B, 2)])
        "Discrete state space system";

    protected
      Integer nx=size(sc.A, 1) "Number of states";
      Integer nu=size(sc.B, 2) "Number of input signals";
      Real LU[nx,nx] "LU decomposition";
      Integer pivots[nx] "Pivots of LU decomposition";
    algorithm
      sd.Ts := Ts;
      sd.method := method;

      if method == Method.ExplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*pre(der_x)
     */
        sd.A := identity(nx) + Ts*sc.A;
        sd.B := Ts*sc.B;
        sd.C := sc.C;
        sd.D := sc.D;
        sd.B2 := zeros(nx, nu);

      elseif method == Method.ImplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*der_x 
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          Ts*sc.A);
        sd.B2 := LU_solve2(
              LU,
              pivots,
              Ts*sc.B);
        sd.A := LU_solve2(
              LU,
              pivots,
              identity(nx));
        sd.B := sd.A*sd.B2;
        sd.C := sc.C;
        sd.D := sd.C*sd.B2 + sc.D;

      elseif method == Method.Trapezoidal then
            /*  der_x = A*x + B*u
             x = pre_x + (Ts/2)*(pre_der_x + der_x); 
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          (Ts/2)*sc.A);
        sd.B2 := LU_solve2(
              LU,
              pivots,
              (Ts/2)*sc.B);
        sd.A := LU_solve2(
              LU,
              pivots,
              identity(nx) + (Ts/2)*sc.A);
        sd.B := sd.A*sd.B2 + sd.B2;
        sd.C := sc.C;
        sd.D := sd.C*sd.B2 + sc.D;

      elseif method == Method.StepExact then
           /* x = phi*pre(x) + gamma*pre(u);
       y = C*x + D*u
    */
        (sd.A,sd.B) := Modelica.Math.Matrices.integralExp(
              sc.A,
              sc.B,
              Ts);
        sd.C := sc.C;
        sd.D := sc.D;
        sd.B2 := zeros(nx, nu);

      elseif method == Method.RampExact then
           /* x = phi*pre(x) + gamma*pre(u) + gamma1/Ts*(u - pre_u);
        -> x = phi*pre(x) + (gamma - gamma1/Ts)*pre(u) + gamma1/Ts*u;
  z = x - gamma1/Ts*u
      leads to
        z = phi*pre(z) + phi*gamma1/Ts*pre(u) +gamma*pre(u) - gamma1/Ts*pre(u) 
       y = C*x + D*u
      -> y = C*z + (D + C*gamma1/Ts)*u
      x = z + gamma1/Ts*u
    
    */
        (sd.A,sd.B,sd.B2) := Modelica.Math.Matrices.integralExpT(
              sc.A,
              sc.B,
              Ts);
        sd.B2 := sd.B2/Ts;
        sd.B := sd.A*sd.B2 + sd.B - sd.B2;
        sd.C := sc.C;
        sd.D := sd.C*sd.B2 + sc.D;

      elseif method == Method.ImpulseExact then
           /* x = phi*pre(x) + phi*B*u;
        y = C*x
        (u = [1,0,1,1,0..,0])
      Limitations: The infinit impulses at t = kT is ignored in the mapping
    */

        sd.A := Modelica.Math.Matrices.exp(sc.A, Ts);
        sd.B := sd.A*sc.B;
        sd.C := sc.C;
        sd.D := sc.C*sc.B;
        sd.B2 := sc.B;

      else
        assert(false, "Argument method (= " + String(method) +
          ") of makeDiscrete is wrong.");
      end if;
      annotation (overloadsConstructor=true);
    end fromStateSpace;

    encapsulated function fromMatrices2
      "Transform a continuous into a discrete linear state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Types.Method;
      import Modelica_LinearSystems2.Math.Matrices.LU_solve2;

      input Real A[:,size(A, 1)] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real B[size(A, 1),:] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real C[:,size(A, 1)] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
    //  input Modelica_LinearSystems2.Types method=Modelica_LinearSystems2.Types.Method.Trapezoidal
      output Modelica_LinearSystems2.DiscreteStateSpace sd(
        redeclare Real A[size(A, 1),size(A, 2)],
        redeclare Real B[size(B, 1),size(B, 2)],
        redeclare Real C[size(C, 1),size(C, 2)],
        redeclare Real D[size(D, 1),size(D, 2)],
        redeclare Real B2[size(B, 1),size(B, 2)]) "Discrete state space system";

    protected
      Integer nx=size(A, 1) "Number of states";
      Integer nu=size(B, 2) "Number of input signals";
      Real LU[nx,nx] "LU decomposition";
      Integer pivots[nx] "Pivots of LU decomposition";
    algorithm
      sd.Ts := Ts;
      sd.method := method;

      if method == Method.ExplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*pre(der_x)
     */
        sd.A := identity(nx) + Ts*A;
        sd.B := Ts*B;
        sd.C := C;
        sd.D := D;
        sd.B2 := zeros(nx, nu);

      elseif method == Method.ImplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*der_x 
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          Ts*A);
        sd.B2 := LU_solve2(
              LU,
              pivots,
              Ts*B);
        sd.A := LU_solve2(
              LU,
              pivots,
              identity(nx));
        sd.B := sd.A*sd.B2;
        sd.C := C;
        sd.D := sd.C*sd.B2 + D;

      elseif method == Method.Trapezoidal then
            /*  der_x = A*x + B*u
             x = pre_x + (Ts/2)*(pre_der_x + der_x); 
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          (Ts/2)*A);
        sd.B2 := LU_solve2(
              LU,
              pivots,
              (Ts/2)*B);
        sd.A := LU_solve2(
              LU,
              pivots,
              identity(nx) + (Ts/2)*A);
        sd.B := sd.A*sd.B2 + sd.B2;
        sd.C := C;
        sd.D := sd.C*sd.B2 + D;

      elseif method == Method.StepExact then
           /* x = phi*pre(x) + gamma*pre(u);
       y = C*x + D*u
    */
        (sd.A,sd.B) := Modelica.Math.Matrices.integralExp(
              A,
              B,
              Ts);
        sd.C := C;
        sd.D := D;
        sd.B2 := zeros(nx, nu);

      elseif method == Method.RampExact then
           /* x = phi*pre(x) + gamma*pre(u) + gamma1/Ts*(u - pre_u);
        -> x = phi*pre(x) + (gamma - gamma1/Ts)*pre(u) + gamma1/Ts*u;
  z = x - gamma1/Ts*u
      leads to
        z = phi*pre(z) + phi*gamma1/Ts*pre(u) +gamma*pre(u) - gamma1/Ts*pre(u) 
       y = C*x + D*u
      -> y = C*z + (D + C*gamma1/Ts)*u
      x = z + gamma1/Ts*u
    
    */
        (sd.A,sd.B,sd.B2) := Modelica.Math.Matrices.integralExpT(
              A,
              B,
              Ts);
        sd.B2 := sd.B2/Ts;
        sd.B := sd.A*sd.B2 + sd.B - sd.B2;
        sd.C := C;
        sd.D := sd.C*sd.B2 + D;

      elseif method == Method.ImpulseExact then
           /* x = phi*pre(x) + phi*B*u;
        y = C*x
        (u = [1,0,1,1,0..,0])
      Limitations: The infinit impulses at t = kT is ignored in the mapping
    */

        sd.A := Modelica.Math.Matrices.exp(A, Ts);
        sd.B := sd.A*B;
        sd.C := C;
        sd.D := C*B;
        sd.B2 := B;

      else
        assert(false, "Argument method (= " + String(method) +
          ") of makeDiscrete is wrong.");
      end if;
      annotation (overloadsConstructor=true);
    end fromMatrices2;
  end 'constructor';

encapsulated operator '-'
    "Contains operators for subtraction of discrete state space systems"

  function subtract
      "Subtraction of two state space systems connected in parallel (= inputs are the same, outputs of the two systems are subtracted)"

      import Modelica;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss1 "State space system 1";
    input DiscreteStateSpace dss2
        "State Space system 2 is subtracted from system 1";
    output DiscreteStateSpace result(
      redeclare Real A[size(dss1.A, 1) + size(dss2.A, 1),size(dss1.A, 2) + size(dss2.A, 2)],
      redeclare Real B[size(dss1.B, 1) + size(dss2.B, 1),size(dss1.B, 2)],
      redeclare Real C[size(dss1.C, 1),size(dss1.C, 2) + size(dss2.C, 2)],
      redeclare Real D[size(dss1.D, 1),size(dss1.D, 2)],
      redeclare Real B2[size(dss1.B2, 1) + size(dss2.B2, 1),size(dss1.B2, 2)])
        "= dss1 - dss2";
    protected
    Integer nx1=size(dss1.A, 1);
    Integer nx2=size(dss2.A, 1);
  algorithm
    assert(abs(dss1.Ts-dss2.Ts)<=Modelica.Constants.eps,"Two discrete state space systems must have the same sample time Ts for subtraction with \"-\".");
    result.A := [dss1.A,zeros(nx1, nx2); zeros(nx2, nx1),dss2.A];
    result.B := [dss1.B; dss2.B];
    result.B2 := [dss1.B2; dss2.B2];
    result.C := [dss1.C,-dss2.C];
    result.D := dss1.D - dss2.D;
    result.Ts := dss1.Ts;
    result.method := dss1.method;
      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> Modelica_LinearSystems2.DiscreteStateSpace.'-'.<b>subtract</b>(dss1, dss2)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This operator function computes the subtraction of two discrete state space systems connected in parallel, i.e. the inputs are the same and the outputs of the two systems are subtracted. Therefore, The systems must have the same number of inputs and outputs but not the same number of states. The resulting system has an order of system_order1 + system_order2.<br>
The operator is used by writing just the following command:
<blockquote><pre>
    dss3 := dss1 - dss2;
</pre></blockquote>

</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   DiscreteStateSpace dss1 = DiscreteStateSpace(A=[-1, 0; 0, -2], B=[1;2], C=[0, 1], D=[0]);
   DiscreteStateSpace dss2 = DiscreteStateSpace(A=[-3, 0; 0, -4], B=[3;4], C=[0, 2], D=[0]);
   
   DiscreteStateSpace dss3;

<b>algorithm</b>
  dss3 := dss1 - dss2;
// dss.A = [-1, 0, 0, 0; 0, -2, 0, 0; 0, 0, -3, 0; 0, 0, 0, -4],
// dss.B = [1; 2; 3; 4],
// dss.C = [0, 1, 0, -2],
// dss.D = [0],
// dss.B2 = [0; 0; 0; 0],
</pre></blockquote>

</html> "));
  end subtract;

  function negate
      "Unary minus (discrete state space system where the output is multiplied by a gain of -1)"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss;
    output DiscreteStateSpace result(
      redeclare Real A[size(dss.A, 1),size(dss.A, 2)],
      redeclare Real B[size(dss.B, 1),size(dss.B, 2)],
      redeclare Real C[size(dss.C, 1),size(dss.C, 2)],
      redeclare Real D[size(dss.D, 1),size(dss.D, 2)]) "= -dss";
  algorithm
    result.A := dss.A;
    result.B := dss.B;
    result.B2 := dss.B2;
    result.C := -dss.C;
    result.D := -dss.D;
    result.Ts := dss.Ts;
    result.method := dss.method;

  end negate;
    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This package contains the <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.'-'.subtract\">'subtract'</a> and the <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.'-'.subtract\">'-'</a> operator for DiscreteStateSpace records. 

</html>"));
end '-';

encapsulated operator function '+'
    "Parallel connection of two discrete state space systems (= inputs are the same, outputs of the two systems are added)"
    import Modelica;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss1 "System 1";
    input DiscreteStateSpace dss2 "System 2 is added in parallel to system 1";
    output DiscreteStateSpace result(
      redeclare Real A[size(dss1.A, 1) + size(dss2.A, 1),size(dss1.A, 2) + size(
        dss2.A, 2)],
      redeclare Real B[size(dss1.B, 1) + size(dss2.B, 1),size(dss1.B, 2)],
      redeclare Real C[size(dss1.C, 1),size(dss1.C, 2) + size(dss2.C, 2)],
      redeclare Real D[size(dss1.D, 1),size(dss1.D, 2)],
      redeclare Real B2[size(dss1.B2, 1) + size(dss2.B2, 1),size(dss1.B2, 2)])
      "= dss1 + dss2";
  protected
    Integer nx1=size(dss1.A, 1);
    Integer nx2=size(dss2.A, 1);
algorithm
    assert(abs(dss1.Ts-dss2.Ts)<=Modelica.Constants.eps,"Two discrete state space systems must have the same sample time Ts for addition with \"+\".");
    result.A := [dss1.A,zeros(nx1, nx2); zeros(nx2, nx1),dss2.A];
    result.B := [dss1.B; dss2.B];
    result.B2 := [dss1.B2; dss2.B2];
    result.C := [dss1.C,dss2.C];
    result.D := dss1.D + dss2.D;
    result.Ts := dss1.Ts;
    result.method := dss1.method;

end '+';

encapsulated operator function '*'
    "Series connection of two discrete state space systems"
    import Modelica;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss1 "System 1";
    input DiscreteStateSpace dss2 "System 2";
    output DiscreteStateSpace result(
      redeclare Real A[size(dss1.A, 1) + size(dss2.A, 1),size(dss1.A, 2) + size(
        dss2.A, 2)],
      redeclare Real B[size(dss1.B, 1) + size(dss2.B, 1),size(dss2.B, 2)],
      redeclare Real C[size(dss1.C, 1),size(dss1.C, 2) + size(dss2.C, 2)],
      redeclare Real D[size(dss1.D, 1),size(dss2.D, 2)],
      redeclare Real B2[size(dss1.B2, 1) + size(dss2.B2, 1),size(dss2.B2, 2)])
      "y = G(s)*u = G(dss1)*G(dss2)*u";
  protected
    Integer nx1=size(dss1.A, 1);
    Integer nx2=size(dss2.A, 1);
algorithm
    if size(dss1.A,1)>0 then
      assert(abs(dss1.Ts-dss2.Ts)<=Modelica.Constants.eps,"Two discrete state space systems must have the same sample time Ts for multiplication with \"*\".");
    end if;
    result.A := [dss1.A,dss1.B*dss2.C; zeros(nx2, nx1),dss2.A];
    result.B := [dss1.B*dss2.D; dss2.B];
    result.B2 := [dss1.B2*dss2.D; dss2.B2];
    result.C := [dss1.C,dss1.D*dss2.C];
    result.D := dss1.D*dss2.D;
    result.Ts := dss2.Ts;
    result.method := dss2.method;

end '*';

encapsulated operator function '=='
    "Check whether two linear discrete state space systems have identical matrices"
    import Modelica;
    import Modelica.Math.Matrices.isEqual;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss1 "System 1";
    input DiscreteStateSpace dss2 "System 2";
    input Real eps(min=0) = 0
      "Two elements e1 and e2 of the two systems are identical if abs(e1-e2) <= eps";
    output Boolean same "=true, if the two systems are identical";
algorithm
    same := isEqual(dss1.A, dss2.A, eps) and isEqual(dss1.B, dss2.B, eps)  and isEqual(dss1.B2, dss2.B2, eps) and isEqual(dss1.C, dss2.C, eps) and isEqual(dss1.D, dss2.D, eps) and abs(dss1.Ts - dss2.Ts) <= Modelica.Constants.eps;
end '==';

// encapsulated operator function 'String'
//     "Transform discrete state space into a String representation"
//    import Modelica;
//    import Modelica_LinearSystems2;
//    import Modelica_LinearSystems2.DiscreteStateSpace;
//    import Modelica.Utilities.Strings;
//
//    input DiscreteStateSpace dss
//       "State space system to be transformed in a String representation";
//    input Integer significantDigits=12
//       "Number of significant digits that are shown";
//    input String name="dss" "Independent variable name used for printing";
//    output String s="";
//
//   protected
//     String space=Strings.repeat(5);
//     String space2=Strings.repeat(3);
//     String space3=Strings.repeat(8);
//     Integer nx=size(dss.A, 1);
//     Integer nu=size(dss.B, 2);
//     Integer ny=size(dss.C, 1);
//     Integer sizeD=size(dss.D, 2);
//     Integer stringMaxLength;
//     Boolean xNamesExist=false;
//     Boolean uNamesExist=false;
//     Boolean yNamesExist=false;
//
// algorithm
//   //Checking if name arrays are empty
//     for i in 1:nx loop
//       xNamesExist := xNamesExist or (dss.xNames[i] <> "");
//     end for;
//
//     for i in 1:ny loop
//       yNamesExist := yNamesExist or (dss.yNames[i] <> "");
//     end for;
//
//     for i in 1:nu loop
//       uNamesExist := uNamesExist or (dss.uNames[i] <> "");
//     end for;
//
//     if xNamesExist then
//       Modelica.Utilities.Streams.print("xNamesExist == true");
//     else
//       Modelica.Utilities.Streams.print("xNamesExist == false");
//     end if;
//
//     stringMaxLength := max(size(dss.xNames, 1), min(size(dss.yNames, 1),
//       11));
//
//     if nx == 0 and sizeD == 0 then
//       s := name + ".A = []\n  " + name + ".B = []\n   " + name + ".C = [] \n   " + name + ".D = []"+ name + ".B2 = []";
//     else
//       s := "\n" + name + ".A = \n";
//
//   //Horizontal
//   // Two alternatives when printing state names
//       if xNamesExist == false then
//         s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
//           "x1 ";
//       else
//         s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
//           dss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(dss.xNames[
//           1]), 11)) + " " + Strings.substring(
//               dss.xNames[1],
//               1,
//               min(Strings.length(dss.xNames[1]), 11));
//       end if;
//
//       for i in 2:nx loop
//
//   //Two alternatives when printing state names
//
//         if xNamesExist == false then
//           s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
//              + String(i - 1))) + "x" + String(i) + " ";
//         else
//           s := s + " " + Strings.repeat(significantDigits + 11 - min(
//             Strings.length(dss.xNames[i - 1]), 11)) + Strings.substring(
//                 dss.xNames[i],
//                 1,
//                 min(Strings.length(dss.xNames[i]), 11));
//
//         end if;
//
//   //s := s + Strings.repeat(6) + "x" + String(i);
//       end for;
//       s := s + "\n";
//
//       for i in 1:nx loop
//   //Vertical
//   //Two alternatives when printing state names
//         if xNamesExist == false then
//           s := s + space + "x" + String(i) + " ";
//         else
//           s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
//             dss.xNames[i]), 11)) + Strings.substring(
//                 dss.xNames[i],
//                 1,
//                 min(Strings.length(dss.xNames[i]), 11)) + " ";
//         end if;
//
//         for j in 1:nx loop
//           if dss.A[i, j] >= 0 then
//             s := s + " ";
//           end if;
//           s := s + String(dss.A[i, j], significantDigits=significantDigits) +
//             Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
//             dss.A[i, j]), significantDigits=significantDigits)));
//        end for;
//         s := s + "\n";
//       end for;
//   //--------------------------------------------------------------------------------------------------------------------------------------------------
//       s := s + "\n" + name + ".B = \n";
//    //Horizontal
//   // Two alternatives when printing state names
//       if uNamesExist == false then
//         s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
//           "u1 ";
//       else
//         s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
//           dss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(dss.uNames[
//           1]), 11)) + " " + Strings.substring(
//               dss.uNames[1],
//               1,
//               min(Strings.length(dss.uNames[1]), 11));
//       end if;
//
//       for i in 2:nu loop
//   //Two alternatives when printing state names
//         if uNamesExist == false then
//           s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
//              + String(i - 1))) + "u" + String(i) + " ";
//         else
//           s := s + " " + Strings.repeat(significantDigits + 11 - min(
//             Strings.length(dss.uNames[i - 1]), 11)) + Strings.substring(
//                 dss.uNames[i],
//                 1,
//                 min(Strings.length(dss.uNames[i]), 11));
//         end if;
//       end for;
//       s := s + "\n";
//   //s := s + Strings.repeat(6) + "x" + String(i);
//       for i in 1:nx loop
//
//   //Vertical
//   //Two alternatives when printing state names
//         if xNamesExist == false then
//           s := s + space + "x" + String(i) + " ";
//         else
//
//           s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
//             dss.xNames[i]), 11)) + Strings.substring(
//                 dss.xNames[i],
//                 1,
//                 min(Strings.length(dss.xNames[i]), 11)) + " ";
//         end if;
//
//         for j in 1:nu loop
//           if dss.B[i, j] >= 0 then
//             s := s + " ";
//           end if;
//           s := s + String(dss.B[i, j], significantDigits=significantDigits) +
//             Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
//             dss.B[i, j]), significantDigits=significantDigits)));
//         end for;
//         s := s + "\n";
//           end for;
//   //--------------------------------------------------------------------------------------------------------------------------------------------------
//       s := s + "\n" + name + ".C = \n";
//    //Horizontal
//   // Two alternatives when printing state names
//       if xNamesExist == false then
//         s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
//           "x1 ";
//       else
//         s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
//           dss.xNames[1]), 11)) + Strings.repeat(min(Strings.length(dss.xNames[
//           1]), 11)) + " " + Strings.substring(
//               dss.xNames[1],
//               1,
//               min(Strings.length(dss.xNames[1]), 11));
//       end if;
//
//       for i in 2:nx loop
//   //Two alternatives when printing state names
//         if xNamesExist == false then
//           s := s + Strings.repeat(significantDigits + 11 - Strings.length("x"
//              + String(i - 1))) + "x" + String(i) + " ";
//         else
//           s := s + " " + Strings.repeat(significantDigits + 11 - min(
//             Strings.length(dss.xNames[i - 1]), 11)) + Strings.substring(
//                 dss.xNames[i],
//                 1,
//                 min(Strings.length(dss.xNames[i]), 11));
//         end if;
//       end for;
//       s := s + "\n";
//   //s := s + Strings.repeat(6) + "x" + String(i);
//
//       for i in 1:ny loop
//   //Vertical
//   //Two alternatives when printing state names
//         if yNamesExist == false then
//           s := s + space + "y" + String(i) + " ";
//         else
//           s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
//             dss.yNames[i]), 11)) + Strings.substring(
//                 dss.yNames[i],
//                 1,
//                 min(Strings.length(dss.yNames[i]), 11)) + " ";
//
//         end if;
//
//         for j in 1:nx loop
//           if dss.C[i, j] >= 0 then
//             s := s + " ";
//           end if;
//           s := s + String(dss.C[i, j], significantDigits=significantDigits) +
//             Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
//             dss.C[i, j]), significantDigits=significantDigits)));
//         end for;
//         s := s + "\n";
//       end for;
//   //--------------------------------------------------------------------------------------------------------------------------------------------------
//       s := s + "\n" + name + ".D = \n";
//    //Horizontal
//   // Two alternatives when printing state names
//       if uNamesExist == false then
//         s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
//           "u1 ";
//       else
//         s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
//           dss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(dss.uNames[
//           1]), 11)) + " " + Strings.substring(
//               dss.uNames[1],
//               1,
//               min(Strings.length(dss.uNames[1]), 11));
//       end if;
//
//       for i in 2:nu loop
//   //Two alternatives when printing state names
//         if uNamesExist == false then
//           s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
//              + String(i - 1))) + "u" + String(i) + " ";
//         else
//           s := s + " " + Strings.repeat(significantDigits + 11 - min(
//             Strings.length(dss.uNames[i - 1]), 11)) + Strings.substring(
//                 dss.uNames[i],
//                 1,
//                 min(Strings.length(dss.uNames[i]), 11));
//         end if;
//       end for;
//       s := s + "\n";
//       for i in 1:ny loop
//   //Vertical
//   //Two alternatives when printing state names
//         if yNamesExist == false then
//           s := s + space + "y" + String(i) + " ";
//         else
//           s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
//             dss.yNames[i]), 11)) + Strings.substring(
//                 dss.yNames[i],
//                 1,
//                 min(Strings.length(dss.yNames[i]), 11)) + " ";
//         end if;
//
//         for j in 1:nu loop
//           if dss.D[i, j] >= 0 then
//             s := s + " ";
//           end if;
//           s := s + String(dss.D[i, j], significantDigits=significantDigits) +
//             Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
//             dss.D[i, j]), significantDigits=significantDigits)));
//         end for;
//         s := s + "\n";
//           end for;
//
// //##################
//
// //--------------------------------------------------------------------------------------------------------------------------------------------------
//       s := s + "\n" + name + ".B2 = \n";
//    //Horizontal
//   // Two alternatives when printing state names
//       if uNamesExist == false then
//         s := s + Strings.repeat(stringMaxLength + significantDigits - 1) +
//           "u1 ";
//       else
//         s := s + Strings.repeat(11 + significantDigits - min(Strings.length(
//           dss.uNames[1]), 11)) + Strings.repeat(min(Strings.length(dss.uNames[
//           1]), 11)) + " " + Strings.substring(
//               dss.uNames[1],
//               1,
//               min(Strings.length(dss.uNames[1]), 11));
//       end if;
//
//       for i in 2:nu loop
//   //Two alternatives when printing state names
//         if uNamesExist == false then
//           s := s + Strings.repeat(significantDigits + 11 - Strings.length("u"
//              + String(i - 1))) + "u" + String(i) + " ";
//         else
//           s := s + " " + Strings.repeat(significantDigits + 11 - min(
//             Strings.length(dss.uNames[i - 1]), 11)) + Strings.substring(
//                 dss.uNames[i],
//                 1,
//                 min(Strings.length(dss.uNames[i]), 11));
//         end if;
//       end for;
//       s := s + "\n";
//   //s := s + Strings.repeat(6) + "x" + String(i);
//       for i in 1:nx loop
//
//   //Vertical
//   //Two alternatives when printing state names
//         if xNamesExist == false then
//           s := s + space + "x" + String(i) + " ";
//         else
//
//           s := s + Strings.repeat(significantDigits + 11 - min(Strings.length(
//             dss.xNames[i]), 11)) + Strings.substring(
//                 dss.xNames[i],
//                 1,
//                 min(Strings.length(dss.xNames[i]), 11)) + " ";
//         end if;
//
//         for j in 1:nu loop
//           if dss.B2[i, j] >= 0 then
//             s := s + " ";
//           end if;
//           s := s + String(dss.B2[i, j], significantDigits=significantDigits) +
//             Strings.repeat(significantDigits + 11 - Strings.length(String(abs(
//             dss.B2[i, j]), significantDigits=significantDigits)));
//         end for;
//         s := s + "\n";
//       end for;
//
// //#################
//
//     end if;
//
//     s := s +"\n\n Ts = " + String(dss.Ts) + "\n method = "+ Modelica_LinearSystems2.Internal.methodString(dss.method);
//
// end 'String';

  encapsulated function timeResponse
    "Compute time response of DiscreteStateSpace system"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace sd "Linear system in discrete state space form";
    input Real u[:,size(sd.B, 2)]=ones(3, size(sd.B, 2))
      "System input (dimension: (input samples) x (number of inputs))";
    input Real x0[size(sd.A, 1)]=zeros(size(sd.A, 1)) "Initial system state";
    output Real y[size(u, 1),size(sd.C, 1)]
      "System response (dimension: (input samples) x (number of outputs))";
    output Real x_continuous[size(u, 1),size(sd.A, 1)]
      "State trajectories (dimension: (input samples) x (number of states)";

  protected
    Integer samples=size(u, 1);
    Integer i;
    Real new_x[size(sd.A, 1)];
    Real x[size(sd.A, 1)]=x0;

  algorithm
    for i in 1:samples loop
      new_x := sd.A*x + sd.B*u[i, :];
      y[i, :] := sd.C*x + sd.D*u[i, :];
      x_continuous[i, :] := x + sd.B2*u[i, :];
      x := new_x;
    end for;
    annotation (Documentation(info="<html>
<p>
Computes the time response of a system in discrete state space form:
</p>
<pre>     <b>x</b>(Ts*(k+1)) = <b>A</b> * <b>x</b>(Ts*k) + <b>B</b> * <b>u</b>(Ts*k)
     <b>y</b>(Ts*k)     = <b>C</b> * <b>x</b>(Ts*k) + <b>D</b> * <b>u</b>(Ts*k)
     <b>x</b>_continuous(Ts*k) = <b>x</b>(Ts*k) + <b>B2</b> * <b>u</b>(Ts*k) 
</pre>
<p>
Note that the system input <b>u</b> must be sampled with the discrete system sample time Ts.
</p>
</html>"));
  end timeResponse;

  encapsulated function initialResponse
    "Compute initial response of DiscreteStateSpace system"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace sd "Linear system in discrete state space form";
    input Real x0[size(sd.A, 1)]=zeros(size(sd.A, 1)) "Initial system state";
    input Integer samples "Number of samples";
    output Real y[samples,size(sd.C, 1)]
      "System response (dimension: (input samples) x (number of outputs))";
    output Real x_continuous[samples,size(sd.A, 1)]
      "State trajectories (dimension: (input samples) x (number of states)";

  protected
    Integer i;
    Real new_x[size(sd.A, 1)];
    Real x[size(sd.A, 1)]=x0;

  algorithm
    for i in 1:samples loop
      new_x := sd.A*x;
      y[i, :] := sd.C*x;
      x_continuous[i, :] := x;
      x := new_x;
    end for;
    annotation (Documentation(info="<html>
<p>
Computes the initial response of a system in discrete state space form:
</p>
<pre>     <b>x</b>(Ts*(k+1)) = <b>A</b> * <b>x</b>(Ts*k)
     <b>y</b>(Ts*k)     = <b>C</b> * <b>x</b>(Ts*k)
     <b>x</b>_continuous(Ts*k) = <b>x</b>(Ts*k) 
</pre>
<p>
Note that the system input <b>u</b> is equal to zero.
</p>
</html>"));
  end initialResponse;

encapsulated package Analysis
    "Functions to analyse discrete state space systems represented by a DiscreteStateSpace record"
encapsulated function eigenValues
      "Calculate the eigenvalues of a linear state space system and write them in a complex vector"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.Math.Complex;

  input DiscreteStateSpace dss "Discrete state space system";
  output Complex eigvalues[size(dss.A, 1)]=Complex.eigenValues(dss.A)
        "eigen values of the system";
algorithm

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  eigenvalues </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>eigenValues</b>(dss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Calculate the eigenvalues of a discrete state space system, i.e. the eigenvalues of the system matrix <b>A</b> of a discrete state space system.
The output is a complex vector containing the eigenvalues.<br>
The eigenvalues <b>ev</b>_d of the discrete system are related to the eigenvalues of the corresponding continuous system <b>ev</b>_c by
</p>
<blockquote>
<p>
<b>ev</b>_d = exp(Ts*<b>ev</b>_c)
</p>
</blockquote>
<p>


</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1;-1,-1],
      B=[1;1],
      C=[1,1],
      D=[0],
      B2=[0;0],
      Ts=0.1);

   DiscreteStateSpace dss=DiscreteStateSpace(ss, Ts=0.1);
 
   Complex eigenvalues[2];
   
<b>algorithm</b>
  eigenvalues = Modelica_LinearSystems2.DiscreteStateSpace.Analysis.eigenValues(dss);
// eigenvalues = {0.900452 + 0.0904977*j, 0.900452 - 0.0904977*j}
//

</pre></blockquote>


</html> "));
end eigenValues;

encapsulated function timeResponse
      "Calculate the time response of a state space system"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Types.TimeResponse;

  input TimeResponse response=TimeResponse.Step;
  extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete(redeclare
          Real y[
             :,size(dss.C, 1),if response == TimeResponse.Initial then 1 else 
      size(dss.B, 2)], redeclare Real x_discrete[:,size(dss.A, 1),if response ==
      TimeResponse.Initial then 1 else size(dss.B, 2)]);// Input/Output declarations of time response functions

  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

    protected
  Real dtVar;
  Real tSpanVar;
  Integer samples;
  Real u[:,size(dss.B, 2)];
//  Real new_x[size(sc.A, 1),1];
//  Real x[size(sc.A, 1),1]=zeros(size(sc.A, 1), 1);

  Real i1;
  Real i2;

algorithm
// set sample time
  if tSpan == 0 then
    tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
  else
    tSpanVar := tSpan;
  end if;

  samples := integer(tSpanVar/dss.Ts + dss.Ts/100 + 1);
// Modelica.Utilities.Streams.print("\nsamples = "+String(samples));

  t := 0:dss.Ts:tSpanVar;
  u := zeros(samples, size(dss.B, 2));

  y := if response == TimeResponse.Initial then zeros(
    samples,
    size(dss.C, 1),
    1) else zeros(
    samples,
    size(dss.C, 1),
    size(dss.B, 2));
  x_discrete := if response == TimeResponse.Initial then zeros(
    samples,
    size(dss.A, 1),
    1) else zeros(
    samples,
    size(dss.A, 1),
    size(dss.B, 2));

  if response == TimeResponse.Initial then
    (y[:, :, 1],x_discrete[:, :, 1]) :=
      Modelica_LinearSystems2.DiscreteStateSpace.Internal.initialResponse1(
      dss,
      x0,
      samples);
  else

    for i1 in 1:size(dss.B, 2) loop
       // Loop over inputs

       // time response to plot
      if response == TimeResponse.Impulse then
        u[1, :] := zeros(size(dss.B, 2));
        u[1, i1] := 1;
      elseif response == TimeResponse.Step then
        u[:, :] := zeros(samples, size(dss.B, 2));
        u[:, i1] := ones(samples);
      elseif response == TimeResponse.Ramp then
        u[:, :] := zeros(samples, size(dss.B, 2));
        u[:, i1] := 0:dss.Ts:tSpanVar + dss.Ts/100;
      else
        assert(false, "Argument response (= " + String(response) + ") of \"Time response to plot\" is wrong.");
      end if;
      (y[:, :, i1],x_discrete[:, :, i1]) :=
        DiscreteStateSpace.Internal.timeResponse1(
        dss,
        u,
        x0);

    end for;
  end if;

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y) </td><td align=center> =  </td>  <td>DiscreteStateSpace.Analysis.<b>timeResponse</b>(responseType, dss)  </td> </tr>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td>DiscreteStateSpace.Analysis.<b>timeResponse</b>(responseType, dss, tSpan, x0)  </td> </tr>

</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function timeResponse calculates the time responses of a discrete state space system. The type of the time response is defined by the input <b>responseType</b>, i.e. 
<blockquote><pre>
    Impulse \"Impulse response\",
    Step \"Step response\",
    Ramp \"Ramp response\",
    Initial \"Initial condition response\"
</pre></blockquote>
Starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and states x are calculated for each time step t=k*dss.Ts.
</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
     A = [0.904837418036],
     B = [0.095162581964],
     C = [2],
     D = [0],
     B2 = [0],
     Ts = 0.1);
  
  Real tSpan= 0.4;
  
  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
  
  Real x0[1]={0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(dss,tSpan,response,x0);
//  y[:,1,1] = {0, 0.1903, 0.3625, 0.5184, 0.6594}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>


</html> "));
end timeResponse;

encapsulated function impulseResponse
      "Calculate the impulse time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
    // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete;
    protected
      Real tSpanVar;
algorithm

// set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=
        Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
        dss=dss,
        tSpan=tSpanVar,
        response=Modelica_LinearSystems2.Types.TimeResponse.Impulse,
        x0=zeros(size(dss.A, 1)));

      annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>impulseResponse</b>(dss)  </td> </tr>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>impulseResponse</b>(dss, tSpan)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and states <b>x</b> are calculated for each time step t=k*dss.Ts.
<blockquote><pre>
DiscreteStateSpace.Analysis.impulseResponse(dss, tSpan)

</pre></blockquote>
gives the same result as
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(dss, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
Note that an appropriate impulse response of a discrete system that is comparable to the impulse response of the corresponding continuous system 
requires the \"ImpulseExact\" conversion from continuous system to discrete system.<br>
See also<br>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a><br>
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.impulseResponse\">StateSpace.Analysis.impulseResponse</a>

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
//  y[:,1,1] = {2, 1.8097, 1.6375, 1.4816, 1.3406}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {1, 0.9048, 0.8187, 0.7408, 0.6703}
</pre></blockquote>


</html> "));
end impulseResponse;

encapsulated function stepResponse
      "Calculate the step time response of a state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
    // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete;
    protected
      Real tSpanVar;
algorithm

// set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) :=
        Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
        dss=dss,
        tSpan=tSpanVar,
        response=Modelica_LinearSystems2.Types.TimeResponse.Step,
        x0=zeros(size(dss.A, 1)));

      annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>stepResponse</b>(dss)  </td> </tr>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>stepResponse</b>(dss, tSpan)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>stepResponse</b> calculates the step response of a discrete state space system. 
Starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and the states <b>x</b> are calculated for each time step t=k*dss.Ts.
<blockquote><pre>
DiscreteStateSpace.Analysis.stepResponse(dss, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(response=Types.TimeResponse.Step, dss, tSpan, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
Note that an appropriate step response of a discrete system that is comparable to the step response of the corresponding continuous system 
requires the \"StepExact\" conversion from continuous system to discrete system.<br>
See also <br>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a><br>
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.stepResponse\">StateSpace.Analysis.stepResponse</a>
</p>


<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A=[0.904837418036 ],
      B=[0.095162581964],
      C=[2],
      D=[0],
      B2=[0],
      Ts=0.1;
      method = Types.Method.StepExact);
      
      

  Real tSpan= 0.4;
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 

<b>algorithm</b>
  (y,t,x) := DiscreteStateSpace.Analysis.stepResponse(dss,tSpan);
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
  import Modelica_LinearSystems2.DiscreteStateSpace;
    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete;
    protected
  Real tSpanVar;
algorithm

// set simulation time span
  if tSpan == 0 then
    tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
  else
    tSpanVar := tSpan;
  end if;

  (y,t,x_discrete) :=
    Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
    dss=dss,
    tSpan=tSpanVar,
    response=Modelica_LinearSystems2.Types.TimeResponse.Ramp,
    x0=zeros(size(dss.A, 1)));

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>rampResponse</b>(dss)  </td> </tr>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>rampResponse</b>(dss, tSpan, x0)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>rampResponse</b> calculates the time response of a discrete state space system for ramp imput u = t. 
Starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dss.Ts.
<blockquote><pre>
DiscreteStateSpace.Analysis.rampResponse(dss, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(response=Types.TimeResponse.Ramp, dss, tSpan, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
Note that an appropriate ramp response of a discrete system that is comparable to the ramp response of the corresponding continuous system 
requires the \"RampExact\" conversion from continuous system to discrete system.<br>

See also<br>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a>
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.rampResponse\">StateSpace.Analysis.rampResponse</a>
</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A=[0.904837418036 ],
      B=[0.095162581964],
      C=[2],
      D=[0.0967483607192],
      B2=[0.0483741803596],
      Ts=0.1;
      method = Types.Method.RampExact);
      
  Real tSpan= 0.4;
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 

<b>algorithm</b>
  (y,t,x) := DiscreteStateSpace.Analysis.rampResponse(dss,tSpan);
//  y[:,1,1] = {0, 0.00967, 0.03746, 0.08164, 0.14064}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0, 0.00484, 0.01873, 0.04082, 0.07032}
</pre></blockquote>

</html> "));
end rampResponse;

encapsulated function initialResponse
      "Calculate the time response of a state space system for given initial condition and zero inputs"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input Real x0[:]=fill(0, 0) "Initial state vector";

    // Input/Output declarations of time response functions:
  extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete(redeclare
          Real y[
             :,size(dss.C, 1),1], redeclare Real x_discrete[:,size(dss.A, 1),1]);
    protected
  Real tSpanVar;

algorithm
  if tSpan == 0 then
    tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
  else
    tSpanVar := tSpan;
  end if;

  (y,t,x_discrete) :=
    Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
    dss=dss,
    tSpan=tSpanVar,
    response=Modelica_LinearSystems2.Types.TimeResponse.Initial,
    x0=x0);

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  (y) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>initialResponse</b>(x0, dss)  </td> </tr>
<tr> <td align=right>  (y, t, x) </td><td align=center> =  </td>  <td> DiscreteStateSpace.Analysis.<b>initialResponse</b>(x0, dss, tSpan)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initialResponse</b> calculates the time response of a discrete state space system for given initial condition and zero inputs. 
tarting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0, the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dss.Ts.
<blockquote><pre>
DiscreteStateSpace.Analysis.initialResponse(x0,dss, dt, tSpan)
</pre></blockquote>
gives the same result as
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(dss, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>
See also <br>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a><br>
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.initialResponse\">StateSpace.Analysis.initialResponse</a><br>
</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A=[0.904837418036 ],
      B=[0.095162581964],
      C=[2],
      D=[0],
      B2=[0],
      Ts=0.1;
      method = Types.Method.StepExact);
      
  Real tSpan= 0.4;
  Real x0[1] = {1};
 
  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1] 


<b>algorithm</b>
  (y,t,x):=DiscreteStateSpace.Analysis.initialResponse(x0,dss,tSpan);
//  y[:,1,1]={2, 1.809, 1.637, 1.4812, 1.3402}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 0.9048, 0.8186, 0.7406, 0.6701}
</pre></blockquote>


</html> "));
end initialResponse;

end Analysis;

encapsulated package Conversion
    "Conversion functions from DiscreteStateSpace into DiscreteTransferFunction"

encapsulated function toDiscreteZerosAndPoles
      "Generate a discrete zeros-and-poles representation from a discrete SISO state space representation"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input DiscreteStateSpace dss "StateSpace object";
  output Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;

    protected
  StateSpace ss=StateSpace(A=dss.A, B=dss.B, C=dss.C, D=dss.D);
  StateSpace ssm=StateSpace.Transformation.toIrreducibleForm(ss);
  Complex poles[:];
  Complex zeros[:];
  Complex cpoles[:];
  Complex czeros[:];
  Real gain;
  Complex frequency;
  Complex cfrequency;
  Complex Gq;
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
    zeros := StateSpace.Internal.invariantZeros2(ssm);
    cpoles := fill(Complex(0),size(poles,1));
    czeros := fill(Complex(0),size(zeros,1));

    if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
      assert(size(ss.B, 2) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
         + String(size(ss.B, 2)) + " instead of 1");
      assert(size(ss.C, 1) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
         + String(size(ss.C, 1)) + " instead of 1");
    end if;
    dzp := DiscreteZerosAndPoles(
        z=zeros,
        p=poles,
        k=1,
        Ts=dss.Ts, method=dss.method);
// set frequency to a complex value which is whether pole nor zero
    for i in 1:size(poles,1) loop
      cpoles[i] := if Complex.'abs'(poles[i])>0 then Complex.log(poles[i])/dss.Ts else Complex(-100);
    end for;
    for i in 1:size(zeros,1) loop
      czeros[i] := if Complex.'abs'(zeros[i])>0 then Complex.log(zeros[i])/dss.Ts else Complex(-100);
    end for;

     v := sum(cat(1, czeros[:].re,  cpoles[:].re))/max(size(czeros,1)+size(cpoles,1),1) + 13/19;
//     v := sum(cat(1, zeros[:].re,  poles[:].re))/max(size(zeros,1)+size(poles,1),1);
    frequency := Complex(v)*17/19;
    cfrequency := Complex.exp(frequency*dss.Ts);
//    cfrequency := frequency;

    Gq := DiscreteZerosAndPoles.Analysis.evaluate(dzp, cfrequency);

    As := -ssm.A;
    for i in 1:size(As, 1) loop
      As[i, i] := As[i, i] + cfrequency.re;
    end for;

    pk := StateSpace.Internal.partialGain(As, ssm.B[:, 1]);
    gain := (ssm.C[1, size(As, 1)]*pk + ss.D[1, 1])/Gq.re;

    dzp := DiscreteZerosAndPoles(
        z=zeros,
        p=poles,
        k=gain,Ts=dss.Ts, method=dss.method);

  else
    dzp := DiscreteZerosAndPoles(
        z=fill(Complex(0), 0),
        p=fill(Complex(0), 0),
        k=scalar(dss.D),Ts=dss.Ts, method=dss.method);

  end if;
//  dzp.uName := dss.uNames[1];
//  dzp.yName := dss.yNames[1];

  annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dzp </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteZerosAndPoles</b>(dss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a DiscreteZerosAndPoles record
 <blockquote><pre>
                 product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
        dzp = k*---------------------------------------------------------
                product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])

</pre></blockquote>of a system from discrete state space representation using the transformation algorithm described in [1].
<br>
The uncontrollable and unobservable parts are isolated and the eigenvalues and invariant zeros of the controllable and observable sub system are calculated.


<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A = [0.9048, 0.0,    0.0;
            0.0,   0.8187, 0.0;
            0.0,   0.0,    0.7408],
      B = [0.09516;         
           0.09063;
           0.0],
      C = [1.0,1.0,1.0],
      D = [0.0],
      Ts = 0.1);
      
 <b>algorithm</b>
  dzp:=Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss);
  
//                         q - 0.860735
//   dzp = 0.1858 -------------------------------
//                 (q - 0.904837)*(q - 0.818731)
</pre></blockquote>


<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V.  </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control,
vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>

</html> "));
end toDiscreteZerosAndPoles;

encapsulated function toDiscreteZerosAndPolesMIMO
      "Generate a zeros-and-poles representation from a MIMO state space representation"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input DiscreteStateSpace dss "DiscreteStateSpace object";

  output DiscreteZerosAndPoles dzp[size(dss.C, 1),size(dss.B, 2)];

    protected
  DiscreteStateSpace dss_siso(
    redeclare Real A[size(dss.A, 1),size(dss.A, 2)],
    redeclare Real B[size(dss.B, 1),1],
    redeclare Real C[1,size(dss.C, 2)],
    redeclare Real D[1,1]);

  Integer ny=size(dss.C, 1);
  Integer nu=size(dss.B, 2);

algorithm
  for ic in 1:ny loop
    for ib in 1:nu loop
      dss_siso := DiscreteStateSpace(
        A=dss.A,
        B=matrix(dss.B[:, ib]),
        C=transpose(matrix(dss.C[ic, :])),
        D=matrix(dss.D[ic, ib]),
        Ts=dss.Ts,
        B2=matrix(dss.B2[:, ib]),
        method=dss.method);
//       dss_siso.uNames := {dss.uNames[ib]};
//       dss_siso.yNames := {dss.yNames[ic]};
//       dss_siso.xNames := dss.xNames;
      dzp[ic, ib] := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(
        dss_siso);
    end for;
  end for;
  annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dzp </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteZerosAndPolesMIMO</b>(dss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a matrix of DiscreteZerosAndPoles records
 <blockquote><pre>
                 product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
       dzp = k*---------------------------------------------------------
                product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre></blockquote>
of a system from discrete state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each single-input-output pair.


<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace ss=Modelica_LinearSystems2.DiscreteStateSpace(
    A = [0.9048, 0.0,    0.0;
          0.0,   0.8187, 0.0;
          0.0,   0.0,    0.7408] 
      
   B = [0.0,      0.09516;         
        0.0906,   0.09063;         
       -0.0864,   0.0],
      
   C = [0.0, 1.0, 1.0;
        1.0, 1.0, 1.0],
      
   D = [1.0, 0.0;
        0.0, 1.0],
        
   B2 = [0.0,  0.0;
         0.0,  0.0;
         0.0,  0.0],
   Ts = 0.1);

<b>algorithm</b>
  dzp:=Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPolesMIMO(dss);

// dzp = [1*(q^2 - 1.55531*q + 0.61012)/( (q - 0.818731)*(q - 0.740818) )
    Ts = 0.1
    method =StepExact,
 0.0906346/(q - 0.818731)
    Ts = 0.1
    method =StepExact;
0.0042407*(q + 0.846461)/( (q - 0.818731)*(q - 0.740818) )
    Ts = 0.1
    method =StepExact,
 1*(q - 0.870319)*(q - 0.667452)/( (q - 0.904837)*(q - 0.818731) )
    Ts = 0.1
    method =StepExact]
    
</pre></blockquote>

<h4><font color=\"#008000\">References</font></h4>
<table>
<tr> <td align=right>  [1] </td><td align=center> Varga, A, Sima, V.  </td>  <td> \"Numerically stable algorithm for transfer function matrix evaluation\"  </td> <td> Int. J. Control,
vol. 33, No. 6, pp. 1123-1133, 1981 </td></tr>
</table>

</html> "));
end toDiscreteZerosAndPolesMIMO;

function toDiscreteTransferFunction
      "Generate a transfer function from a SISO state space representation"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteTransferFunction;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  input DiscreteStateSpace dss "DiscreteStateSpace object";

  output DiscreteTransferFunction dtf "DiscreteTransferFunction object";

    protected
  DiscreteZerosAndPoles dzp=DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss);

algorithm
  dtf := DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp);

  annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dtf </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteTransferFunction</b>(dss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a DiscreteTransferFunction record
<blockquote><pre>
           n(z)     b0 + b1*z + ... + bn*z^n
  dtf = -------- = -------------------------- 
           d(z)     a0 + a1*z + ... + an*z^n
 </pre></blockquote>

The algorithm uses <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles\">toDiscreteZerosAndPoles</a> to convert the
discrete state space system into a discrete zeros and poles representation first and after that <a href=\"Modelica://Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction\">ZerosAndPoles.Conversion.toTransferFunction</a> to generate the transfer function.



<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A = [0.9048, 0.0,    0.0;
            0.0,   0.8187, 0.0;
            0.0,   0.0,    0.7408],
      B = [0.09516;         
           0.09063;
           0.0],
      C = [1.0,1.0,1.0],
      D = [0.0],
      B2 = [0, 0;
            0, 0;
            0, 0],
      Ts = 0.1);
      
 <b>algorithm</b>
  dtf:=Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss);
  
//             0.185797*z - 0.159922
//   dtf = -------------------------------
//           z^2 - 1.72357*z + 0.740818
</pre></blockquote>



</html> "));
end toDiscreteTransferFunction;

function toDiscreteTransferFunctionMIMO
      "Generate a discrete transfer function of a MIMO system from discrete state space representation"
  import Modelica_LinearSystems2;

  import Modelica;
  import Modelica_LinearSystems2.DiscreteTransferFunction;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  input DiscreteStateSpace dss "DiscreteStateSpace object";

  output DiscreteTransferFunction dtf[size(dss.C, 1),size(dss.B, 2)]
        "Matrix of discrete transfer function objects";

    protected
  DiscreteZerosAndPoles dzp[:,:];
  parameter Integer m=size(dss.B, 2);
  parameter Integer p=size(dss.C, 1);

algorithm
  dzp := Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPolesMIMO(dss);
  for i1 in 1:m loop
    for i2 in 1:p loop
      dtf[i2, i1] := DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp[i2, i1]);
    end for;
  end for;

      annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dtf </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteTransferFunctionMIMO</b>(dss)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Computes a matrix of DiscreteTransferFunction records
<blockquote><pre>
           n(z)     b0 + b1*z + ... + bn*z^n
  dtf = -------- = -------------------------- 
           d(z)     a0 + a1*z + ... + an*z^n
 </pre></blockquote>
with repetitive application of <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteTransferFunction\">Conversion.toDiscreteTransferFunction</a>


<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
   A = [0.9048, 0.0,    0.0;
          0.0,   0.8187, 0.0;
          0.0,   0.0,    0.7408] 
      
   B = [0.0,      0.09516;         
        0.0906,   0.09063;         
       -0.0864,   0.0],
      
   C = [0.0, 1.0, 1.0;
        1.0, 1.0, 1.0],
      
   D = [1.0, 0.0;
        0.0, 1.0],
        
   B2 = [0.0,  0.0;
         0.0,  0.0;
         0.0,  0.0],
   Ts = 0.1);

<b>algorithm</b>
  dtf:=Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss);

// dtf = [(1*z^2 - 1.55531*z + 0.61012)/(z^2 - 1.55955*z + 0.606531)
 
 Ts = 0.1
 method =StepExact,
 0.0906346/(z - 0.818731)
 
 Ts = 0.1
 method =StepExact;
(0.0042407*z + 0.00358958)/(z^2 - 1.55955*z + 0.606531)
 
 Ts = 0.1
 method =StepExact,
 (1*z^2 - 1.53777*z + 0.580896)/(z^2 - 1.72357*z + 0.740818)
 
 Ts = 0.1
 method =StepExact]
</pre></blockquote>

</html> "));
end toDiscreteTransferFunctionMIMO;

end Conversion;

encapsulated package Design
function kfStepMatrices
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
extends Modelica.Icons.Function;

import Modelica;
import Modelica_LinearSystems2;
import Modelica_LinearSystems2.Math;
import Modelica_LinearSystems2.Math.Matrices.LAPACK;
import Modelica_LinearSystems2.DiscreteStateSpace;

input Real A[:,size(A, 1)] "Transition matrix of the discrete system";
input Real B[size(A, 1),:] "Input matrix of the discrete system";
input Real C[:,size(A, 1)] "Output matrix of the discrete system";
input Real P[size(A, 1),size(A, 1)]
        "State covariance matrix of the previous instant";
input Real Q[size(B, 2),size(B, 2)]
        "Input or process noise covariance matrix of the previous instant";
input Real R[size(C, 1),size(C, 1)]
        "Output or measurement noise covariance matrix of the previous instant";

output Real K[size(A, 1),size(C, 1)]=P*transpose(C) "Kalman filter gain matrix";
output Real P_new[size(A, 1),size(A, 1)] "Updated state covariance matrix";
output Real UMutri[size(C, 1),size(C, 1)]
        "Square root (left Cholesky factor) of the covariance matrix M";
output Real M[size(C, 1),size(C, 1)]=DiscreteStateSpace.Internal.symMatMul(C, P, R, true)
        "Upper triangle of measurement prediction covariance C*P*C' + R";

output Real PCT[size(A, 1),size(C, 1)]=P*transpose(C) "Matrix P*C'";
    protected
Integer nx=size(A, 1) "Number of states, order of the system";
Integer nu=size(B, 2) "Number of inputs";
Integer ny=size(C, 1) "number of outputs";
Integer l1;
Integer l2;
Real alpha=1.0;

// Real P[size(A, 1),size(A, 1)]=P_in    "State covariance matrix of the previous instant";

Integer info;

algorithm
//PCT:=P*transpose(C) "Matrix P*C'";
//M:=DiscreteStateSpace.Internal.symMatMul(C, P, R, true)     "Upper triangle of measurement prediction covariance C*P*C' + R";
(UMutri,info) := LAPACK.dpotrf(M, true);// Calculate the Cholesky factorization U*U' of M
assert(info == 0, "Calculating a Cholesky decomposition with function \"Lapack.dpotrf\" failed in function \"kfStep\".");
for l1 in 2:ny loop
  for l2 in 1:l1 - 1 loop
    UMutri[l1, l2] := 0.0;
  end for;
end for;

      //K from K*M = P*C' with K*U'*U = P*C', U is Cholesky factor
K := LAPACK.dtrsm(UMutri, K, alpha, true, true, false, false);
K := LAPACK.dtrsm(UMutri, K, alpha, true, true, true, false);

// Calculate upper triangle of symmetric P-K*C*P
for l1 in 1:nx loop
  for l2 in l1:nx loop
    P_new[l1, l2] := P[l1, l2] - K[l1, :]*PCT[l2, :];
  end for;
end for;
//Calculate upper triangle of A*(P-K*C*P)*A'
P_new := DiscreteStateSpace.Internal.symMatMul(A, P_new, P_new, false);
//Calculate upper triangle of A*(P-K*C*P)*A' + B*Q*B'
P_new := DiscreteStateSpace.Internal.symMatMul(B, Q, P_new, true);

      // Note that P_new contains the upper triangle of the symmetric covariance matrix.
      // To complete the matrix, the strict lower triangle could be calculated by
for l1 in 2:nx loop
  for l2 in 1:l1-1 loop
    P_new[l1,l2] := P_new[l2,l1];
  end for;
end for;

end kfStepMatrices;

function kfStepState
      "One step, i.e.estimation of the state vector using a kalman filter iteration for discrete systems"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input DiscreteStateSpace dss;
  input Real P[size(dss.A, 1),size(dss.A, 1)]
        "State covariance matrix of the previous instant";
  input Real Q[size(dss.B, 2),size(dss.B, 2)]
        "Input or process noise covariance matrix of the previous instant";
  input Real R[size(dss.C, 1),size(dss.C, 1)]
        "Output or measurement noise covariance matrix of the previous instant";
  input Real x[size(dss.A, 1)] "Estimated state vector of previous instant";
  input Real u[size(dss.B, 2)] "input vector";
  input Real y[size(dss.C, 1)] "Measured output vector";

  output Real x_new[size(dss.A, 1)];
  output Real K[size(dss.A, 1),size(dss.C, 1)] "Kalman filter gain matrix";
  output Real P_new[size(dss.A, 1),size(dss.A, 1)]
        "Updated state covariance matrix";

    protected
  Real UMutri[size(dss.C, 1),size(dss.C, 1)]
        "Square root (left Cholesky factor) of the covariance matrix M=R + C*P*C'";
   Real z[size(dss.A, 1)];

algorithm
  (K,P_new,UMutri) := DiscreteStateSpace.Design.kfStepMatrices(
    dss.A,
    dss.B,
    dss.C,
    P,
    Q,
    R);
  z := dss.A*x + dss.B*u;
  x_new := z - K*(dss.C*z - y);

end kfStepState;

end Design;

encapsulated package Plot
    "Functions to plot discrete state space system responses"
encapsulated function bodeSISO
      "Plot bode plot of the corresponding discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  input DiscreteStateSpace dss "discrete state space system";
  input Integer iu=1 "index of input";
  input Integer iy=1 "index of output";
  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min=0.1
        "Minimum frequency value, if autoRange = false";
  input Modelica.SIunits.Frequency f_max=10
        "Maximum frequency value, if autoRange = false";

  input Boolean magnitude=true "= true, to plot the magnitude of dtf" 
                                                                     annotation(choices(__Dymola_checkBox=true));
  input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(__Dymola_checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

    protected
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp
        "Zeros and poles Transfer functions to be plotted";
  DiscreteStateSpace dss_siso(
    redeclare Real A[size(dss.A, 1),size(dss.A, 2)],
    redeclare Real B[size(dss.B, 1),1],
    redeclare Real C[1,size(dss.C, 2)],
    redeclare Real D[1,1],
    redeclare Real B2[size(dss.B2, 1),1]);

algorithm
  assert(iu <= size(dss.B, 2) and iu > 0, "index for input is " + String(iu) + " which is not in [1, "
     + String(size(dss.B, 2)) + "].");
  assert(iy <= size(dss.C, 1) and iy > 0, "index for output is " + String(iy) + " which is not in [1, "
     + String(size(dss.C, 1)) + "].");
  dss_siso := DiscreteStateSpace(
    A=dss.A,
    B=matrix(dss.B[:, iu]),
    C=transpose(matrix(dss.C[iy, :])),
    D=matrix(dss.D[iy, iu]),
    B2=matrix(dss.B2[:, iu]),
    Ts=dss.Ts,
    method=dss.method);
  dzp := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss_siso);

  DiscreteZerosAndPoles.Plot.bode(
    dzp,
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
DiscreteStateSpace.Plot.<b>plotBodeSISO</b>(dss)
   or
DiscreteStateSpace.Plot.<b>plotBodeSISO</b>(dss, iu, iy, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
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
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
   A = [0.9048, 0.0,    0.0;
          0.0,   0.8187, 0.0;
          0.0,   0.0,    0.7408] 
      
   B = [0.0;         
        0.0906;         
       -0.0864],
      
   C = [0.0, 1.0, 1.0],
      
   D = [1.0],
        
   B2 = [0.0;
         0.0;
         0.0],         
   Ts = 0.1);

   
   Integer iu=1;
   Integer iy=1;


<b>algorithm</b>
   Modelica_LinearSystems2.DiscreteStateSpace.Plot.plotBodeSISO(dss, iu, iy)
//  gives:
</pre></blockquote>

</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagDis.png\">
</p>
<p>
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhaseDis.png\">
</p>
<p>


</html> "));
end bodeSISO;

encapsulated function timeResponse
      "Plot the time response of the system. The response type is selectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteStateSpace dss;
  input Real tSpan=0 "Simulation time span [s]";

  input TimeResponse response=TimeResponse.Step;

  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

  input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Time response"));

    protected
  Plot.Records.Curve curve;
  Integer i1;
  Integer i2;
  Plot.Records.Diagram diagram2[size(dss.C, 1)];

  Real y[:,size(dss.C, 1),if response == TimeResponse.Initial then 1 else size(dss.B,2)]
        "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  Real t[:] "Time vector: (number of samples)";
  Real x[:,size(dss.A, 1),if response == TimeResponse.Initial then 1 else size(dss.B,2)]
        "State trajectories: (number of samples) x (number of states) x (number of inputs)";

  Real yy[:,:,:] "Output response";
  Real tt[:] "Time vector: (number of samples)";

  String yNames[size(dss.C, 1)];
  String uNames[size(dss.B, 2)];
  Integer loops=if response == TimeResponse.Initial then 1 else size(dss.B,2);

algorithm
  (y,t,x) := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
    dss=dss,
    tSpan=tSpan,
    response=response,
    x0=x0);

  tt := fill(0,2*size(t,1)-1);
  yy := fill(0,2*size(t,1)-1,size(y,2),size(y,3));

  for i in 1:size(t,1)-1 loop
    tt[2*i-1] := t[i];
    tt[2*i] := t[i+1];
    yy[2*i-1,:,:] := y[i,:,:];
    yy[2*i,:,:] := y[i,:,:];
  end for;
  tt[size(tt,1)] := t[size(t,1)];
  yy[size(tt,1),:,:] := y[size(t,1),:,:];

// generate headings
  for i1 in 1:size(dss.B, 2) loop
    uNames[i1] := "u" + String(i1);
  end for;
  for i1 in 1:size(dss.C, 1) loop
    yNames[i1] := "y" + String(i1);
  end for;

  // for i1 in 1:size(dss.B, 2) loop
  //   uNames[i1] := if dss.uNames[i1] == "" then "u" + String(i1) else dss.uNames[
  //     i1];
  // end for;
  // for i1 in 1:size(dss.C, 1) loop
  //   yNames[i1] := if dss.yNames[i1] == "" then "y" + String(i1) else dss.yNames[
  //     i1];
  // end for;

  for i2 in 1:loops loop
    for i1 in 1:size(dss.C, 1) loop
      curve := Plot.Records.Curve(
        x=tt,
        y=yy[:, i1, i2],
        autoLine=true);

      diagram2[i1] := defaultDiagram;
      diagram2[i1].curve := {curve};
      diagram2[i1].heading := if response == TimeResponse.Initial then defaultDiagram.heading +" "+ yNames[i1] else defaultDiagram.heading + "  " + uNames[i2] + " -> " + yNames[i1];
      diagram2[i1].yLabel := yNames[i1];

    end for;

    if subPlots then
      Plot.diagramVector(diagram2, device);
    else
      for i1 in 1:size(dss.C, 1) loop
        Plot.diagram(diagram2[i1], device);
      end for;
    end if;
  end for;

  annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>timeResponse</b>(dss);
or
DiscreteStateSpace.Plot.<b>timeResponse</b>(dss, tSpan,response, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>


<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>timeResponse</b> plots the time response of a discrete state space system. The character of the time response if defined by the input <tt>response</tt>, i.e. Impulse, Step, Ramp, or Initial. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.plotImpulse\">plotImpulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[0.0,1.0,1.0],
D=[0.0])

Real Ts = 0.1;
Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
DiscreteStateSpace dss=DiscreteStateSpace(ss,Ts,method);
Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

<b>algorithm</b>
Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse(dss, response=response);

</pre></blockquote>
</html> "));
end timeResponse;

encapsulated function impulse "Impulse response plot"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";
    input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

    input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Impulse response"));

    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Impulse
        "type of time response";
    Real tSpanVar;

algorithm
      // set sample time
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
    else
      tSpanVar := tSpan;
    end if;

    Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse(
      dss=dss,
      tSpan=tSpanVar,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>impulse</b>(dss);
or
DiscreteStateSpace.Plot.<b>impulse</b>(dss, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>plotImpulse</b> plots the impulse responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[0.0,1.0,1.0],
D=[0.0])

Real Ts = 0.1;
Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.ImpulseExact;
DiscreteStateSpace dss=DiscreteStateSpace(dss,Ts,method);

<b>algorithm</b>
Modelica_LinearSystems2.StateSpace.Plot.impulse(dss)

</p>
</pre></blockquote>

</html> "));
end impulse;

encapsulated function step "Step response plot"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";

    input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Step response"));

    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Step "type of time response";

    Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

    Real tSpanVar;

algorithm
// set sample time
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
    else
      tSpanVar := tSpan;
    end if;

    Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse(
      dss=dss,
      tSpan=tSpanVar,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>step</b>(dss);
or
DiscreteStateSpace.Plot.<b>step</b>(dss, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>step</b> plots the discrete step responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initial\">initial</a>.





</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[0.0,1.0,1.0],
D=[0.0])

Real Ts = 0.1;
Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
DiscreteStateSpace dss=DiscreteStateSpace(dss,Ts,method);

<b>algorithm</b>
Modelica_LinearSystems2.StateSpace.Plot.step(dss, tSpan=3)
</pre></blockquote>

</html> "));
end step;

encapsulated function ramp "Ramp response plot"

    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.Types.TimeResponse;

    import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";

    input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Ramp response"));
    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Ramp "type of time response";

    Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

    Real tSpanVar;

algorithm
// set sample time
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
    else
      tSpanVar := tSpan;
    end if;

    DiscreteStateSpace.Plot.timeResponse(
      dss=dss,
      tSpan=tSpanVar,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>ramp</b>(ss);
or
DiscreteStateSpace.Plot.<b>ramp</b>(dss, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>ramp</b> plots the ramp responses of a discrete state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initial\">initial</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[1.0,1.0,1.0],
D=[0.0])

Real Ts = 0.1;
Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.RampExact;
DiscreteStateSpace dss=DiscreteStateSpace(dss,Ts,method);

<b>algorithm</b>
Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp(dss)
</pre></blockquote>


</html> "));
end ramp;

encapsulated function initialResponse "Initial condition response plot"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.Types.TimeResponse;

  import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteStateSpace dss;
  input Real tSpan=0 "Simulation time span [s]";
  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

  input Boolean subPlots=true
        "true if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(Dialog,choices(__Dymola_checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Initial response"));
    protected
  input Modelica_LinearSystems2.Types.TimeResponse response=
      Modelica_LinearSystems2.Types.TimeResponse.Initial
        "type of time response";

    Real tSpanVar;

algorithm
// set sample time
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
    else
      tSpanVar := tSpan;
    end if;
  Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse(
    dss=dss,
    tSpan=tSpanVar,
    response=response,
    x0=x0,
    defaultDiagram=defaultDiagram,
    device=device);

  annotation (interactive=true, Documentation(info="<html> 
<p><b><font style=\"color: #008000; \">Syntax</font></b></p>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>initial</b>(ss);
or
DiscreteStateSpace.Plot.<b>initial</b>(dss, tSpan, x0, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Function <b>initial</b> plots the initial responses of a discrete state space system for the initial state vector x0 for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>, 
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>, and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>.



</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
A=[-1.0,0.0,0.0; 0.0,-2.0,3.0; 0.0,-2.0,-3.0],
B=[1.0; 1.0; 0.0],
C=[0.0,1.0,1.0],
D=[0.0])

Real x0={1,0.5,0.5}; 

Real Ts = 0.1;
Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
DiscreteStateSpace dss=DiscreteStateSpace(dss,Ts,method);


<b>algorithm</b>
Modelica_LinearSystems2.DiscreteStateSpace.Plot.initial(dss, x0=x0)
</pre></blockquote>
</html> "));
end initialResponse;

end Plot;

encapsulated package Import

function fromModel
      "Generate a StateSpace data record by linearization of a model"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.DiscreteStateSpace;

  input String modelName "Name of the Modelica model" annotation(Dialog(translatedModel));
  input Real T_linearize=0 "point in time of simulation to linearize the model";
  input String fileName="dslin" "Name of the result file";
  input Modelica.SIunits.Time Ts=1 "Sample time";
  input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";

    protected
  String fileName2=fileName + ".mat";
  Boolean OK1=simulateModel(
      problem=modelName,
      startTime=0,
      stopTime=T_linearize);
  Boolean OK2=importInitial("dsfinal.txt");
  Boolean OK3=linearizeModel(
      problem=modelName,
      resultFile=fileName,
      startTime=T_linearize,
      stopTime=T_linearize + 3*Ts);

  Real nxMat[1,1]=readMatrix(
      fileName2,
      "nx",
      1,
      1);
  Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
  Integer nx=integer(nxMat[1, 1]);
  Integer nu=ABCDsizes[2] - nx;
  Integer ny=ABCDsizes[1] - nx;
  Real ABCD[nx + ny,nx + nu]=readMatrix(
      fileName2,
      "ABCD",
      nx + ny,
      nx + nu);
  String xuyName[nx + nu + ny]=readStringMatrix(
      fileName2,
      "xuyName",
      nx + nu + ny);

  StateSpace ss(
    redeclare Real A[nx,nx],
    redeclare Real B[nx,nu],
    redeclare Real C[ny,nx],
    redeclare Real D[ny,nu]) "= model linearized at initial point";
    public
  output DiscreteStateSpace result(
    redeclare Real A[nx,nx],
    redeclare Real B[nx,nu],
    redeclare Real B2[nx,nu],
    redeclare Real C[ny,nx],
    redeclare Real D[ny,nu]) "= discrete model linearized at initial point";

algorithm
  ss.A := ABCD[1:nx, 1:nx];
  ss.B := ABCD[1:nx, nx + 1:nx + nu];
  ss.C := ABCD[nx + 1:nx + ny, 1:nx];
  ss.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
  ss.uNames := xuyName[nx + 1:nx + nu];
  ss.yNames := xuyName[nx + nu + 1:nx + nu + ny];
  ss.xNames := xuyName[1:nx];

  result := DiscreteStateSpace(
    sc=ss,
    Ts=Ts,
    method=method);

  annotation (interactive=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> DiscreteStateSpace.Import.<b>fromModel</b>(modelName, T_linearize, fileName, Ts, method)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Generate a DiscreteStateSpace data record by linearization of a model defined by modelName. The linearization is performed at time T_linearize of the simulation.
The result of linearization is transformed into a StateSpace record and then converted into a DiscreteStateSpace record

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   String modelName = \"Modelica_LinearSystems2.Examples.Utilities.DoublePendulum\"; 
   Real T_linearize = 5;

<b>algorithm</b>
  dss = Modelica_LinearSystems2.DiscreteStateSpace.Import.fromModel(modelName, T_linearize);

// ss.A  = [1, 0.1, 0.0338415578929376, 0.00169207789464688, -0.010114371760331, -0.000505718588016548;
            0, 1, 0.676831157858752, 0.0338415578929376, -0.202287435206619,   -0.010114371760331;
            0, 0, 0.892698457060726, 0.0946349228530364, -0.0895698633812754,  -0.00447849316906376;
            0, 0, -2.14603085878547, 0.892698457060726, -1.79139726762551, -0.0895698633812755;
            0, 0, 0.0738077919657481, 0.0036903895982874, 1.0110083777752, 0.10055041888876;
            0, 0, 1.47615583931496, 0.0738077919657481, 0.220167555503916, 1.0110083777752];

// ss.B= [0.00170323086692055;
          0.0165800882263292;
         -0.00215003506298196;
         -0.02103518146498;
          0.00152042385523347;
          0.0144812915324601];
          
// ss.C=identity(6),
// ss.D= [0.000437113227802044;
          0.00874226455604088;
         -0.000549137994799829;
         -0.0109827598973295;
          0.000398179639293453;
          0.00796359278610463];
          
ss.B2  = [0.000437113227802044;
          0.00874226455604088;
         -0.000549137994866478;
         -0.0109827598973295;
          0.000398179639305232;
          0.00796359278610463];
      

                
</pre></blockquote>



</html> 
"));
end fromModel;

  encapsulated function fromFile "Read a StateSpace data record from mat-file"

    import Modelica;
    import Modelica_LinearSystems2.DiscreteStateSpace;
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
    output DiscreteStateSpace result(
      redeclare Real A[nx,nx],
      redeclare Real B[nx,nu],
      redeclare Real B2[nx,nu],
      redeclare Real C[ny,nx],
      redeclare Real D[ny,nu]) "= model linearized at initial point";

    protected
    Real ABCD[nx + ny,nx + nu]=Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
        fileName,
        matrixName,
        nx + ny,
        nx + nu);
    Real B2[nx,nu]=Modelica_LinearSystems2.Internal.Streams.readMatrixInternal(
        fileName,
        "B2",
        nx,
        nu);
    Real Ts[1,1]=readMatrix(
        fileName,
        "Ts",
        1,
        1);

  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.B := B2;
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    result.Ts := scalar(Ts);
    Modelica.Utilities.Streams.print("StateSpace record loaded from file: \"" +
      fileName + "\"");

    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> DiscreteStateSpace.Import.<b>fromFile</b>(fileName, matrixName)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Reads and loads a discrete state space system from a mat-file <tt>fileName</tt>. The file must contain the matrix [A, B; C, D] named matrixName, the matris B2 and the integer nx representing the order of the system, i.e. the number of rows of the square matrix A.

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
     

<b>algorithm</b>
  dss:=Modelica_LinearSystems2.DiscreteStateSpace.Import.fromFile(\"dss.mat\");
//  dss=StateSpace(
      A=[-4.5, 1.5, 4.0; -4.0, 1.0, 4.0; -1.5, -0.5, 1],
      B=[2; 1; 2],
      C=[1, 0, 0],
      D=[0],
      B2=[0;0;0],
      Ts=0.2,
      method=Trapezoidal);


</pre></blockquote>


</html> "));
  end fromFile;

end Import;

encapsulated package Internal
function symMatMul
      "Calculate the upper triangle of A*B*A'+a*C with B and C symmetric"
  extends Modelica.Icons.Function;
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,:];
  input Real B[size(A, 2),size(A, 2)];
  input Real C[size(A, 1),size(A, 1)];
  input Boolean add=true "true if a==1, false if a==0";
  output Real M[size(A, 1),size(A, 1)];

    protected
  Integer a1=size(A, 1);
  Integer a2=size(A, 2);
  Integer i;
  Integer j;

  Real alpha=1.0;
  Real beta=if add then 1 else 0;

  Real Butri[a2,a2]=B;
  Real Cutri[a1,a1]=C;
  Real Ah[a1,a2];

algorithm
  for i in 1:a2 loop
    Butri[i, i] := B[i, i]/2;
  end for;
  if add then
    for i in 1:a1 loop
      Cutri[i, i] := C[i, i]/2;
    end for;
    for i in 2:a1 loop
      for j in 1:i - 1 loop
        Cutri[i, j] := 0.0;
      end for;
    end for;
  end if;

  Ah := LAPACK.dtrmm(Butri, A, alpha, true, true, false, false);
  M := LAPACK.dgemm(Ah, A, Cutri, alpha, beta, false, true);
// M:= Ah*transpose(A)+Cutri;
  for i in 1:a1 loop
    for j in i:a1 loop
      M[i,j] := M[i,j]+M[j,i];
    end for;
  end for;

end symMatMul;

function timeResponseSamples
      "Estimate reasonable discretisation sample time and simulation time span for time response plot"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Complex;

  input Modelica_LinearSystems2.DiscreteStateSpace dss;
  output Real tSpan "Time span";
    protected
  Modelica_LinearSystems2.Math.Complex eig[size(dss.A, 1)];
  Real realp[size(dss.A, 1)];
  Real sorted[size(dss.A, 1)];
  Real indices[size(dss.A, 1)];
  Integer i;
algorithm
  eig := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.eigenValues(dss);
  for i in 1:size(dss.A, 1) loop
    eig[i] := if Complex.'abs'(eig[i])>1e-10 then Complex.log(eig[i])/dss.Ts else Complex(-100);
  end for;

  //eig := Complex.log(eig)/dss.Ts;
  for i in 1:size(eig, 1) loop
    realp[i] := eig[i].re;
  end for;
  (sorted,indices) := Modelica.Math.Vectors.sort(realp);

  // Estimate simulation time span
  if sorted[end] < 0 then
    tSpan := -5/sorted[end];
  elseif sorted[end] > 0 then
    tSpan := 15/sorted[end];
  elseif sorted[end] == 0 then
    tSpan := 15000;
  end if;
  // Curb simulation time span to maximal 15000s
  if tSpan > 15000 then
    tSpan := 15000;
  end if;

end timeResponseSamples;

  encapsulated function initialResponse1
      "Compute initial response of DiscreteStateSpace system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss "Linear system in discrete state space form";
    input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial system state";
    input Integer samples "Number of samples";
    output Real y[samples,size(dss.C, 1)]
        "System response (dimension: (input samples) x (number of outputs))";
    output Real x[samples,size(dss.A, 1)]
        "State trajectories (dimension: (input samples) x (number of states)";

    protected
    Integer i;
    Real new_x[size(dss.A, 1)];
    Real xi[size(dss.A, 1)]=x0;

  algorithm
    for i in 1:samples loop
      new_x := dss.A*xi;
      y[i, :] := dss.C*xi;
      x[i, :] := xi;
      xi := new_x;
    end for;
    annotation (Documentation(info="<html>
<p>
Computes the initial response of a system in discrete state space form:
</p>
<pre>     <b>x</b>(Ts*(k+1)) = <b>A</b> * <b>x</b>(Ts*k)
     <b>y</b>(Ts*k)     = <b>C</b> * <b>x</b>(Ts*k)
     <b>x</b>_continuous(Ts*k) = <b>x</b>(Ts*k) 
</pre>
<p>
Note that the system input <b>u</b> is equal to zero.
</p>
</html>"));
  end initialResponse1;

encapsulated function timeResponse1
      "Compute time response of DiscreteStateSpace system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

  input DiscreteStateSpace dss "Linear system in discrete state space form";
  input Real u[:,size(dss.B, 2)]=ones(3, size(dss.B, 2))
        "System input (dimension: (input samples) x (number of inputs))";
  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial system state";
  output Real y[size(u, 1),size(dss.C, 1)]
        "System response (dimension: (input samples) x (number of outputs))";
  output Real x[size(u, 1),size(dss.A, 1)]
        "State trajectories (dimension: (input samples) x (number of states)";

    protected
  Integer samples=size(u, 1);
  Integer i;
  Real new_x[size(dss.A, 1)];
  Real xi[size(dss.A, 1)]=x0;

algorithm
  for i in 1:samples loop
    new_x := dss.A*xi + dss.B*u[i, :];
    y[i, :] := dss.C*xi + dss.D*u[i, :];
    x[i, :] := xi;
    xi := new_x;
  end for;
  annotation (Documentation(info="<html>
<p>
Computes the time response of a system in discrete state space form:
</p>
<pre>     <b>x</b>(Ts*(k+1)) = <b>A</b> * <b>x</b>(Ts*k) + <b>B</b> * <b>u</b>(Ts*k)
     <b>y</b>(Ts*k)     = <b>C</b> * <b>x</b>(Ts*k) + <b>D</b> * <b>u</b>(Ts*k)
     <b>x</b>_continuous(Ts*k) = <b>x</b>(Ts*k) + <b>B2</b> * <b>u</b>(Ts*k) 
</pre>
<p>
Note that the system input <b>u</b> must be sampled with the discrete system sample time Ts.
</p>
</html>"));
end timeResponse1;
end Internal;

  annotation (
    defaultComponentName="stateSpaceDiscrete",
    Documentation(info="<html>
<p>
This record defines a linear time invariant difference
equation system in state space form:
</p>
<pre>     <b>x</b>(Ts*(k+1)) = <b>A</b> * <b>x</b>(Ts*k) + <b>B</b> * <b>u</b>(Ts*k)
     <b>y</b>(Ts*k)     = <b>C</b> * <b>x</b>(Ts*k) + <b>D</b> * <b>u</b>(Ts*k)
     <b>x</b>_continuous(Ts*k) = <b>x</b>(Ts*k) + <b>B2</b> * <b>u</b>(Ts*k) 
</pre>
<p>
with
</p>
<ul>
<li> <b>Ts</b> - the sample time</li>
<li> <b>k</b> - the index of the actual sample instance (k=0,1,2,3,...)</li>
<li> <b>t</b> - the time</li>
<li> <b>u</b>(t) - the input vector,</li>
<li> <b>y</b>(t) - the output vector,</li>
<li> <b>x</b>(t) - the discrete state vector (x(t=Ts*0) is the initial state),</li>
<li> <b>x</b>_continuous(t) - the state vector of the continuous system
     from which the discrete block has been derived (details see below),</li>
<li> <b>A,B,C,D,B2</b> - matrices of appropriate dimensions.</li>
</ul>
<p>
A discrete system is usually derived by discretization from a 
continuous block, e.g., by function
LinearSystems.DiscreteStateSpace.fromStateSpace.
If the discretization method, e.g., the trapezoidal method,
accesses <b>actual and past</b> values of the input <b>u</b> 
(e.g. <b>u</b>(Ts*k), <b>u</b>(Ts*(k-1), <b>u</b>(Ts*(k-2))), 
a state transformation is needed to get the difference equation 
above where only the actual value <b>u</b>(Ts*k) is accessed. 
</p>
<p>
If the original continuous state vector should be computed
from the sampled data system above, the matrices of this 
transformation have to be known. For simplicity and efficiency,
here only the specific transformation used by function
LinearSystems.DiscreteStateSpace.fromStateSpace is stored in 
the data record of the discrete system via matrix <b>B2</b>.
Therefore, the state vector of the underlying continuous
system can be calculated by adding the term <b>B2</b>*<b>u</b> to the
state vector of the discretized system.
</p>
<p>
In Modelica notation, the difference equation above
can be implemented as:
</p>
<pre>
     <b>when</b> {<b>initial</b>(), <b>sample</b>(Ts,Ts)} <b>then</b>
        new_x = A * x + B * u;
            y = C * x + D * u;
            x = <b>pre</b>(new_x);
        x_continuous = x + B2 * u;
     <b>end when</b>;
</pre>
<p>
Since no \"next value\" operator is available in Modelica, an
auxiliary variable new_x stores the value of x for the
next sampling instant. The relationship between new_x and x
is defined via equation \"x = <b>pre</b>(new_x)\".
</p>
<p>
The body of the when-clause is active during initialization and at the
next sample instant t=Ts. Note, the when-equation is not
active after the initialization at t=0 (due to <b>sample</b>(Ts,Ts)), 
since the state x of the initialization has to be used also at t=0.
</p>
<p>
In library Blocks.<b>Controller</b> additional equations are
added for the initialization to uniquely compute the
initial vector x:
</p>
<pre>
  <b>initial equation</b> 
     <b>if</b> init == InitialState <b>then</b>
        x = x_start;
     <b>elseif</b> init == SteadyState <b>then</b>
        x = new_x;
     <b>end if</b>;
</pre>
<p>
Optionally, x is set to a given start vector x_start.
As <b>default initialization</b>, the equation \"x = new_x\" is
added that defines steady state initialization for the
discrete system. As a consequence, the output y(Ts*k), k=0,1,2,..,
remains constant after this initialization, 
provided the input vector u(Ts*k) remains constant.
</p>
</html>"));
end DiscreteStateSpace;
