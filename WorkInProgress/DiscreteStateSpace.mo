within Modelica_LinearSystems2.WorkInProgress;
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

  encapsulated operator 'constructor'
    import Modelica_LinearSystems2;
  function fromDiscreteTransferFunction = 
      Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteStateSpace
                                                                                                      annotation (Documentation(info="<html> </html>"));
  function fromDiscreteZerosAndPoles = 
      Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteStateSpace
                                                                                                      annotation (Documentation(info="<html> </html>"));
    function fromMatrices "Default constructor for a DiscreteStateSpace record"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
      output Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace sd(
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

      function fromReal "Generate a StateSpace data record from a Real value"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
      output Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace sd(
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

  encapsulated function timeResponse
    "Compute time response of DiscreteStateSpace system"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
    import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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

encapsulated package Conversion
    "Conversion functions from DiscreteStateSpace into DiscreteTransferFunction"
function toDiscreteTransferFunction
      "Generate a transfer function from a SISO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

  input DiscreteStateSpace dss "DiscreteStateSpace object";

  output Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "DiscreteTransferFunction object";

    protected
  StateSpace ss = StateSpace(dss.A,dss.B,dss.C,dss.D);
  TransferFunction tf = StateSpace.Conversion.toTransferFunction(ss);

algorithm
  dtf := DiscreteTransferFunction(n=tf.n, d=tf.d, Ts=dss.Ts, method=dss.method);

    annotation (Documentation(info="<html>




</html> "));
end toDiscreteTransferFunction;

encapsulated function toDiscreteZerosAndPoles
      "Generate a discrete zeros-and-poles representation from a discrete SISO state space representation"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
      cpoles[i] := Complex.log(poles[i])/dss.Ts;
    end for;
    for i in 1:size(zeros,1) loop
      czeros[i] := Complex.log(zeros[i])/dss.Ts;
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
end toDiscreteZerosAndPoles;

end Conversion;

encapsulated package Plot
    "Functions to plot discrete state space system responses"
encapsulated function bodeSISO
      "Plot bode plot of the corresponding discrete transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

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
  Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "Transfer functions to be plotted";
  DiscreteStateSpace dss_siso(
    redeclare Real A[size(dss.A, 1),size(dss.A, 2)],
    redeclare Real B[size(dss.B, 1),1],
    redeclare Real C[1,size(dss.C, 2)],
    redeclare Real D[1,1],
    redeclare Real B2[size(dss.B, 1),1]);

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
  dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss_siso);

  DiscreteTransferFunction.Plot.bode(
    dtf,
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
StateSpace.Plot.<b>plotBodeSISO</b>(ss, iu, iy, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
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

end Plot;

encapsulated package Internal
function kfStepMatrices
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
Real mnorm;
Real mrcond;
Real eps=1e-12;

// Real P[size(A, 1),size(A, 1)]=P_in    "State covariance matrix of the previous instant";

Integer info;

algorithm
//PCT:=P*transpose(C) "Matrix P*C'";
//M:=DiscreteStateSpace.Internal.symMatMul(C, P, R, true)     "Upper triangle of measurement prediction covariance C*P*C' + R";
mnorm := LAPACK.dlansy(M, "1", true);
(UMutri,info) := LAPACK.dpotrf(M, true);// Calculate the Cholesky factorization U*U' of M
assert(info == 0, "Calculating a Cholesky decomposition with function \"Lapack.dpotrf\" failed in function \"kfMatrices\".");
(mrcond, info) := LAPACK.dpocon(UMutri, mnorm, true);
assert(info == 0, "Calculating the reciprocal condition number with function \"Lapack.dpocon\" failed in function \"kfMatrices\".");
assert(mrcond>eps, "Matrix C*P*C' + R in function \"kfMatrices\" is (numerically) singular");

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
P_new := symmetric(P_new);

end kfStepMatrices;

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

function kfStepState
      "One step, i.e.estimation of the state vector using a kalman filter iteration for discrete systems"
  extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
  (K,P_new,UMutri) := DiscreteStateSpace.Internal.kfStepMatrices(
    dss.A,
    dss.B,
    dss.C,
    P,
    Q,
    R);
  z := dss.A*x + dss.B*u;
  x_new := z - K*(dss.C*z - y);

end kfStepState;

function kfStepMatrices2
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;

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
output Real R_new[ny,ny];
    protected
Integer nx=size(A, 1) "Number of states, order of the system";
Integer nu=size(B, 2) "Number of inputs";
Integer ny=size(C, 1) "number of outputs";
Integer l1;
Integer l2;
Real alpha=1.0;

// Real P[size(A, 1),size(A, 1)]=P_in    "State covariance matrix of the previous instant";

Real ICK[ny,ny];
Real Iy[ny,ny]=identity(ny);

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

ICK:= (Iy-C*K)*transpose(Iy-C*K);
R_new := Modelica.Math.Matrices.solve2(ICK*R,Iy-ICK);
M:=DiscreteStateSpace.Internal.symMatMul(
        C,
        P,
        R,
        true);
(UMutri,info) := LAPACK.dpotrf(M, true);
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

end kfStepMatrices2;

end Internal;

encapsulated package Design "Design functions for DiscreteStateSpace"

function sr_kfStepMatrices
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
  extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,size(A, 1)] "Transition matrix of the discrete system";
  input Real C[:,size(A, 1)] "Output matrix of the discrete system";
  input Real S[size(A, 1),size(A, 1)]
        "Cholesky factor of the state covariance matrix of the previous instant";
  input Real wB[size(A, 1),:]
        "Weighting matrix of the input noise covariance matrix of the previous instant";
  input Real Cq[size(wB,2),size(wB,2)]
        "Cholesky factor of input or process noise covariance matrix of the previous instant";
  input Real Cr[size(C, 1),size(C, 1)]
        "Cholesky factor or measurement noise covariance matrix of the previous instant";

  output Real K[size(A, 1),size(C, 1)] "Kalman filter gain matrix";
  output Real S_new[size(A, 1),size(A, 1)] "Updated state covariance matrix";
  output Real Cr_new[size(C, 1),size(C, 1)]
        "Modified Cholesky output or measurement noise covariance matrix";
    output Real M[size(C, 1)+size(A, 1), 2*size(C, 1)+size(A, 1)];

    protected
  Integer nx=size(A,1);
  Integer ny=size(C,1);

  Integer info;

algorithm
  (M,,info) := LAPACK.dgelqf([Cr, C*S, zeros(ny,ny); zeros(nx,ny), A*S, wB*Cq]);
  assert(info == 0, "Computation of LQ factorization with \"Lapack.dgelqf\" failed in function \"sr_kfStepMatrices\".");
  S_new := Matrices.triangle(M[ny+1:ny+nx,ny+1:nx+ny],false);
  Cr_new := Matrices.triangle(M[1:ny,1:ny],false);
  K := Matrices.Internal.solve2rSym(Cr_new,M[ny+1:ny+nx,1:ny],true,false);
end sr_kfStepMatrices;

function sr_kfStepMatrices2
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
  extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,size(A, 1)] "Transition matrix of the discrete system";
  input Real C[:,size(A, 1)] "Output matrix of the discrete system";
  input Real S[size(A, 1),size(A, 1)]
        "Cholesky factor of the state covariance matrix of the previous instant";
//  input Real Cq[size(A,1),size(A,1)] "Cholesky factor of input or process noise covariance matrix of the previous instant";
  input Real Cq[size(A, 1),size(A, 1)]
        "Cholesky factor of input or process noise covariance matrix of the previous instant";
  input Real Cr[size(C, 1),size(C, 1)]
        "Cholesky factor or measurement noise covariance matrix of the previous instant";

  output Real K[size(A, 1),size(C, 1)] "Kalman filter gain matrix";
  output Real S_new[size(A, 1),size(A, 1)] "Updated state covariance matrix";
  output Real Cr_new[size(C, 1),size(C, 1)]
        "Modified Cholesky output or measurement noise covariance matrix";
//    output Real M[size(C, 1)+size(A, 1), size(C, 1)+2*size(A, 1)];
    output Real M[size(C, 1)+size(A, 1), size(C, 1)+size(A, 1) +size(Cq,2)];

    protected
  Integer nx=size(A,1);
  Integer ny=size(C,1);
  Integer nu=size(Cq,2);

  Integer info;

algorithm
//  (M,,info) := LAPACK.dgelqf([Cr, C*S, zeros(ny,nx); zeros(nx,ny), A*S, Cq]);
//M:=[Cr, C*S, zeros(ny,nu); zeros(nx,ny), A*S, Cq];
  (M,,info) := LAPACK.dgelqf([Cr, C*S, zeros(ny,nu); zeros(nx,ny), A*S, Cq]);
  assert(info == 0, "Computation of LQ factorization with \"Lapack.dgelqf\" failed in function \"sr_kfStepMatrices\".");
  for i in 1:nx+ny loop
    if M[i,i]<0 then
      M[:,i]:=-M[:,i];
    end if;
  end for;

  S_new := Matrices.triangle(M[ny+1:ny+nx,ny+1:nx+ny],false);
  Cr_new := Matrices.triangle(M[1:ny,1:ny],false);
  K := Matrices.Internal.solve2rSym(Cr_new,M[ny+1:ny+nx,1:ny],true,false);
//  K := M[ny+1:ny+nx,1:ny];
end sr_kfStepMatrices2;

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

  output Real K[size(A, 1),size(C, 1)]=P*transpose(C)
        "Kalman filter gain matrix";
  output Real P_new[size(A, 1),size(A, 1)] "Updated state covariance matrix";

   output Real UMutri[size(C, 1),size(C, 1)]
        "Square root (left Cholesky factor) of the covariance matrix M";
  output Real M[size(C, 1),size(C, 1)];//=DiscreteStateSpace.Internal.symMatMul(C, P, R, true) "Upper triangle of measurement prediction covariance C*P*C' + R";

  output Real PCT[size(A, 1),size(C, 1)];//=P*transpose(C) "Matrix P*C'";

    protected
  Integer nx=size(A, 1) "Number of states, order of the system";
  Integer nu=size(B, 2) "Number of inputs";
  Integer ny=size(C, 1) "number of outputs";
  Integer l1;
  Integer l2;
  Real alpha=1.0;
  Real mnorm;
  Real mrcond;

  Real eps=1e-5;

// Real P[size(A, 1),size(A, 1)]=P_in    "State covariance matrix of the previous instant";

  Integer info;

algorithm
PCT:=P*transpose(C) "Matrix P*C'";
M:=Math.Matrices.Internal.symMatMul(     C, P, R, true)
        "Upper triangle of measurement prediction covariance C*P*C' + R";
  mnorm := LAPACK.dlansy(M, "1", true);
  (UMutri,info) := LAPACK.dpotrf(M, true);// Calculate the Cholesky factorization U*U' of M
  assert(info == 0, "Calculating a Cholesky decomposition with function \"Lapack.dpotrf\" failed in function \"kfMatrices\".");
  (mrcond,info) := LAPACK.dpocon(UMutri, mnorm, true);
  assert(info == 0, "Calculating the reciprocal condition number with function \"Lapack.dpocon\" failed in function \"kfMatrices\".");
//  assert(mrcond > eps, "Matrix C*P*C' + R in function \"kfMatrices\" is (numerically) singular");

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
  P_new := Math.Matrices.Internal.symMatMul(     A, P_new, P_new, false);
//Calculate upper triangle of A*(P-K*C*P)*A' + B*Q*B'
  P_new := Math.Matrices.Internal.symMatMul(     B, Q, P_new, true);

  // Note that P_new contains the upper triangle of the symmetric covariance matrix.
  // To complete the matrix, the strict lower triangle could be calculated by
  P_new := symmetric(P_new);

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
  (K,P_new,UMutri) := Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Design.kfStepMatrices(
                                                               dss.A, dss.B, dss.C, P, Q, R);
  z := dss.A*x + dss.B*u;
  x_new := z - K*(dss.C*z + dss.D*u - y);

end kfStepState;

function kfStepMatrices2
      "One step, i.e. prediction and update of a kalman filter iteration for discrete systems"
  extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;
      import Modelica_LinearSystems2.DiscreteStateSpace;

  input Real A[:,size(A, 1)] "Transition matrix of the discrete system";
//  input Real B[size(A, 1),:] "Input matrix of the discrete system";
  input Real C[:,size(A, 1)] "Output matrix of the discrete system";
  input Real P[size(A, 1),size(A, 1)]
        "State covariance matrix of the previous instant";
  input Real Q[size(A, 1),size(A, 1)]
        "Input or process noise covariance matrix of the previous instant";
  input Real R[size(C, 1),size(C, 1)]
        "Output or measurement noise covariance matrix of the previous instant";

  output Real K[size(A, 1),size(C, 1)] "Kalman filter gain matrix";
  output Real P_new[size(A, 1),size(A, 1)] "Updated state covariance matrix";

    protected
  Real M[size(C, 1),size(C, 1)];//=DiscreteStateSpace.Internal.symMatMul(C, P, R, true) "Upper triangle of measurement prediction covariance C*P*C' + R";

  Real PCT[size(A, 1),size(C, 1)];//=P*transpose(C) "Matrix P*C'";

  Integer nx=size(A, 1) "Number of states, order of the system";
//  Integer nu=size(B, 2) "Number of inputs";
  Integer ny=size(C, 1) "number of outputs";
  Integer l1;
  Integer l2;
  Real alpha=1.0;

  Integer info;

algorithm
  PCT:=P*transpose(C) "Matrix P*C'";
  M:=Math.Matrices.Internal.symMatMul_C(     C, P, R, true);
//  M:=C*PCT+R;

//   mnorm := LAPACK.dlansy(M, "1", true);
//   (UMutri,info) := LAPACK.dpotrf(M, true);// Calculate the Cholesky factorization U*U' of M
//   assert(info == 0, "Calculating a Cholesky decomposition with function \"Lapack.dpotrf\" failed in function \"kfMatrices\".");
//   (mrcond,info) := LAPACK.dpocon(UMutri, mnorm, true);
//   assert(info == 0, "Calculating the reciprocal condition number with function \"Lapack.dpocon\" failed in function \"kfMatrices\".");
//   assert(mrcond > eps, "Matrix C*P*C' + R in function \"kfMatrices\" is (numerically) singular");

  K := Math.Matrices.Internal.solve2rSym_C(M, PCT,false, true);
//  K := Math.Matrices.solve2rSym(M, PCT,false, true);

// Calculate upper triangle of symmetric P-K*C*P
  for l1 in 1:nx loop
    for l2 in l1:nx loop
      P_new[l1, l2] := P[l1, l2] - K[l1, :]*PCT[l2, :];
    end for;
  end for;

////Calculate upper triangle of A*(P-K*C*P)*A'
//  P_new := DiscreteStateSpace.Internal.symMatMul_C(A, P_new, P_new, false);
//  P_new := P_new + Q;
//Calculate upper triangle of A*(P-K*C*P)*A' + Q
//   P_new := symmetric(P_new);
//   P_new := A*P_new*transpose(A)+Q;
  P_new := Math.Matrices.Internal.symMatMul_C(A, P_new, Q, true);

  // Note that P_new contains the upper triangle of the symmetric covariance matrix.
  // To complete the matrix, the strict lower triangle could be calculated by
  P_new := symmetric(P_new);

end kfStepMatrices2;

function sr_ukfEstimate_2
      "Calculate filter gain and the updated mean of the state and Cholesky factor S of covariance P of a UKF"
  extends Modelica.Icons.Function;

      import Modelica;
      import Modelica_LinearSystems2;
//  import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Math.Matrices;
      import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real y[:] "Measured output vector";
  input Real xm[:] "Predicted state mean";
  input Real ym[size(y,1)] "Predicted mean";
  input Real S[size(xm,1),size(xm, 1)] "Cholesky factor of covariance matrix";
  input Real Syy[size(y,1),size(y, 1)]
        "Cholesky factor of transformed covariance matrix";
  input Real Pxy[size(xm, 1),size(y, 1)] "Transformed cross covariance matrix";

  output Real K[size(xm,1),size(y,1)] "Filter gain";
  output Real Su[size(xm, 1),size(xm, 1)]
        "Updated Cholesky factor state covariance matrix";
  output Real xmu[size(xm, 1)] "Updated state mean";

    protected
  Real U[size(xm,1),size(y,1)];
algorithm
  K := LAPACK.dtrsm(Syy, Pxy, 1, true, false, true, false);
  K := LAPACK.dtrsm(Syy, K, 1, true, false, false, false);

//  U := K*transpose(Syy);
  U := K*Syy;
  Su := Matrices.choleskyDownDate(S,U[:,1],false);
   for i in 2:size(y,1) loop
     Su := Matrices.choleskyDownDate(Su,U[:,i],false);
   end for;

//  K[2,:]:={0.01,0.01};
  xmu := xm + K*(y-ym);

end sr_ukfEstimate_2;
end Design;
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
