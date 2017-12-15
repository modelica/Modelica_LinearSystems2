within Modelica_LinearSystems2;
operator record DiscreteStateSpace
  "Discrete state space description of a linear, time invariant difference equation system (data + operations)"

  Real A[:,size(A, 1)]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real B[size(A, 1),:]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real C[:,size(A, 1)]  annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));
  Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="new_x = A*x + B*u;  y = C*x + D*u;  x_cont = x + B2*u"));

  Modelica.SIunits.Time Ts=1 "Sample time"
       annotation(Dialog(group="Data used to construct discrete from continuous system"));
  Real B2[size(B, 1),size(B, 2)]=fill(0,size(B,1),size(B,2))
    "Reconstruct continuous state"
       annotation(Dialog(group="Data used to construct discrete from continuous system"));
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));

//      String yNames[size(C, 1)]=fill("", size(C, 1)) "Names of the output signals" annotation(Dialog(group="Signal names"));
//      String xNames[size(A, 1)]=fill("", size(A, 1)) "Names of the states"  annotation(Dialog(group="Signal names"));
//      String uNames[size(B, 2)]=fill("", size(B, 2)) "Names of the input signals" annotation(Dialog(group="Signal names"));

  encapsulated operator 'constructor'
    "Collection of operators to construct a DiscreteStateSpace data record"
    import Modelica_LinearSystems2;
    import Modelica;

    function fromDiscreteTransferFunction =
      Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteStateSpace
      "Generate a DiscreteStateSpace data record from a discrete transfer function"
    annotation (Documentation(info="<html> </html>"));

    function fromDiscreteZerosAndPoles =
      Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteStateSpace
      "Generate a DiscreteStateSpace data record from a discrete zeros-and-poles description"
    annotation (Documentation(info="<html> </html>"));

    encapsulated function fromReal
      "Generate a DiscreteStateSpace data record from a real value"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
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
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = 'constructor'.<b>fromReal</b>(r)
   or
dss = 'constructor'.<b>fromReal</b>(r, Ts, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function constructs a DiscreteStateSpace record dss from a Real value,
i.e. a discrete state space system without a state and an output without dynamics:
</p>
<blockquote><pre>
y = r*u
</pre></blockquote>
<p>
Therefore, the matrices are defined by
</p>
<blockquote><pre>
dss.A = fill(0,0,0);
dss.B = fill(0,0,1);
dss.C = fill(0,1,0);
dss.D = [r];
dss.B2 = fill(0,0,1);
</pre></blockquote>
<p>
The default values of sample time <b>Ts</b> and discretization method <b>method</b> are
</p>
<blockquote><pre>
    Ts = 1
method = Modelica_LinearSystems2.Types.Method.Trapezoidal
</pre></blockquote>
<p>
respectively.
</p>
</html>"));
    end fromReal;

    function fromMatrices
      "Generate a DiscreteStateSpace data record from A, B, C and D matrices"
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
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization methodDiscretization method" annotation (Dialog(group="Data used to construct discrete from continuous system"));
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
<h4>Syntax</h4>
<blockquote>
<pre>
dss = 'constructor'.<b>fromMatrices</b>(A, B, C, D)
dss = 'constructor'.<b>fromMatrices</b>(A, B, C, D, Ts, B2, method)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This function constructs a DiscreteStateSpace record dss with
</p>
<blockquote><pre>
dss.A = A;
dss.B = B;
dss.C = C;
dss.D = D;
dss.B2 = B2;
dss.Ts = Ts;
dss.method = method;
</pre></blockquote>
<p>
i.e. the input-matrices are the system matrices of the discrete system.
The default values of sample time <b>Ts</b> and discretization method
<b>method</b> are
</p>
<blockquote><pre>
    Ts = 1
method = Modelica_LinearSystems2.Types.Method.Trapezoidal
</pre></blockquote>
<p>
respectively. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.'constructor'.fromMatrices2\">fromMatrices2</a>
where the inputs are the matrices of a continuous system which is to convert to discrete state space.
</p>

<h4>Example</h4>
<blockquote><pre>
  import dss=Modelica_LinearSystems2.DiscreteStateSpace;
  Real A[1,1] = [1];
  Real B[1,1] = [1];
  Real C[1,1] = [1];
  Real D[1,1] = [0];

public
  DiscreteStateSpace dss;

<b>algorithm</b>
  dss := dss.'constructor'.fromMatrices(A, B, C, D)  // or just: dss := dss(A, B, C, D, Ts=1, B2=[0], method=Modelica_LinearSystems2.Types.Method.Trapezoidal);
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
      "Generate a DiscreteStateSpace data record from a continuous state space system "
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Utilities.Types.Method;
      import Modelica_LinearSystems2.Math.Matrices.LU_solve2;

      input Modelica_LinearSystems2.StateSpace ss
        "Continuous linear state space system";
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
      output Modelica_LinearSystems2.DiscreteStateSpace dss(
        redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
        redeclare Real B[size(ss.B, 1),size(ss.B, 2)],
        redeclare Real C[size(ss.C, 1),size(ss.C, 2)],
        redeclare Real D[size(ss.D, 1),size(ss.D, 2)],
        redeclare Real B2[size(ss.B, 1),size(ss.B, 2)])
        "Discrete state space system";

    protected
      Integer nx=size(ss.A, 1) "Number of states";
      Integer nu=size(ss.B, 2) "Number of input signals";
      Real LU[nx,nx] "LU decomposition";
      Integer pivots[nx] "Pivots of LU decomposition";
    algorithm
      dss.Ts := Ts;
      dss.method := method;

      if method == Method.ExplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*pre(der_x)
     */
        dss.A := identity(nx) + Ts*ss.A;
        dss.B := Ts*ss.B;
        dss.C := ss.C;
        dss.D := ss.D;
        dss.B2 := zeros(nx, nu);

      elseif method == Method.ImplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*der_x
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          Ts*ss.A);
        dss.B2 := LU_solve2(
              LU,
              pivots,
              Ts*ss.B);
        dss.A := LU_solve2(
              LU,
              pivots,
              identity(nx));
        dss.B := dss.A*dss.B2;
        dss.C := ss.C;
        dss.D := dss.C*dss.B2 + ss.D;

      elseif method == Method.Trapezoidal then
            /*  der_x = A*x + B*u
             x = pre_x + (Ts/2)*(pre_der_x + der_x);
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          (Ts/2)*ss.A);
        dss.B2 := LU_solve2(
              LU,
              pivots,
              (Ts/2)*ss.B);
        dss.A := LU_solve2(
              LU,
              pivots,
              identity(nx) + (Ts/2)*ss.A);
        dss.B := dss.A*dss.B2 + dss.B2;
        dss.C := ss.C;
        dss.D := dss.C*dss.B2 + ss.D;

      elseif method == Method.StepExact then
           /* x = phi*pre(x) + gamma*pre(u);
       y = C*x + D*u
    */
        (dss.A,dss.B) := Modelica.Math.Matrices.integralExp(
              ss.A,
              ss.B,
              Ts);
        dss.C := ss.C;
        dss.D := ss.D;
        dss.B2 := zeros(nx, nu);

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
        (dss.A,dss.B,dss.B2) := Modelica.Math.Matrices.integralExpT(
              ss.A,
              ss.B,
              Ts);
        dss.B2 := dss.B2/Ts;
        dss.B := dss.A*dss.B2 + dss.B - dss.B2;
        dss.C := ss.C;
        dss.D := dss.C*dss.B2 + ss.D;

      elseif method == Method.ImpulseExact then
           /* x = phi*pre(x) + phi*B*u;
        y = C*x
        (u = [1,0,1,1,0..,0])
      Limitations: The infinite impulses at t = kT is ignored in the mapping
    */

        dss.A := Modelica.Math.Matrices.exp(ss.A, Ts);
        dss.B := dss.A*ss.B;
        dss.C := ss.C;
        dss.D := ss.C*ss.B;
        dss.B2 := ss.B;

      else
        assert(false, "Argument method (= " + String(method) +
          ") of makeDiscrete is wrong.");
      end if;
      annotation (
        Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = 'constructor'.<b>fromStateSpace</b>(ss, Ts)
dss = 'constructor'.<b>fromStateSpace</b>(ss, Ts, method)
</pre></blockquote>

<h4>Description</h4>
<p>
This function derives a linear time invariant difference
equation system in state space form
</p>
<blockquote><pre>
<b>x</b>(Ts*(k+1))        = <b>A</b> * <b>x</b>(Ts*k) + <b>B</b> * <b>u</b>(Ts*k)
<b>y</b>(Ts*k)            = <b>C</b> * <b>x</b>(Ts*k) + <b>D</b> * <b>u</b>(Ts*k)
<b>x</b>_continuous(Ts*k) =     <b>x</b>(Ts*k) + <b>B2</b> * <b>u</b>(Ts*k)
</pre></blockquote>
<p>
with
</p>
<ul>
<li> <b>Ts</b> - the sample time,</li>
<li> <b>k</b> - the index of the actual sample instance (k=0,1,2,3,...),</li>
<li> <b>t</b> - the time,</li>
<li> <b>u</b>(t) - the input vector,</li>
<li> <b>y</b>(t) - the output vector,</li>
<li> <b>x</b>(t) - the discrete state vector (x(t=Ts*0) is the initial state),</li>
<li> <b>x</b>_continuous(t) - the state vector of the continuous system
     from which the discrete block has been derived (details see below),</li>
<li> <b>A, B, C, D, B2</b> - matrices of appropriate dimensions.</li>
</ul>
<p>
from continuous state space form
</p>
<blockquote><pre>
der(<b>xc</b>(t)) = <b>ss.A</b> * <b>xc</b>(t) + <b>ss.B</b> * <b>us</b>(t)
    <b>yc</b>(t)  = <b>ss.C</b> * <b>xc</b>(t) + <b>ss.D</b> * <b>uc</b>(t)
</pre></blockquote>
<p>
The applied discretization method is selected by the user from
</p>
<ul>
<li> <b>ExplicitEuler</b> - Discretization with explicit Euler integration,</li>
<li> <b>ImplicitEuler</b> - Discretization with implicit Euler integration,</li>
<li> <b>Trapezoidal</b> - Discretization with trapezoidal integration (Tustins method, recommended),</li>
<li> <b>ImpulseExact</b> - Exact discretization for impulse inputs,</li>
<li> <b>StepExact</b> - Exact discretization for step inputs (zero-order hold equivalent),</li>
<li> <b>RampExact</b> - Exact discretization for ramp inputs (first-order hold equivalent).</li>
</ul>

<h4>Example</h4>
<blockquote><pre>
  import dss=Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.StateSpace;

  StateSpace ss=StateSpace(A = [1], B = [1], C = [1], D = [0]);
  Modelica.SIunits.Time Ts=0.1;
  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal;

public
  DiscreteStateSpace dss;

<b>algorithm</b>
  dss := dss.'constructor'.fromStateSpace(ss, Ts);

  //or just:
  //dss := dss(ss=ss, Ts=Ts, method=method);

  //  dss.A = [1.1053],
  //  dss.B = [0.11080],
  //  dss.C = [1],
  //  dss.D = [0.0526],
  //  dss.Ts = 0.1,
  //  dss.B2 = [0.0526],
  //  dss.method = Modelica_LinearSystems2.Types.Method.Trapezoidal
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.'constructor'.fromMatrices2\">fromMatrices2</a>
</p>
</html>"));
    end fromStateSpace;

    encapsulated function fromMatrices2
      "Generate a DiscreteStateSpace data record from matrices of a continuous state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Utilities.Types.Method;
      import Modelica_LinearSystems2.Math.Matrices.LU_solve2;

      input Real A[:,size(A, 1)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
      input Real B[size(A, 1),:] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
      input Real C[:,size(A, 1)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
      input Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
      input Modelica.SIunits.Time Ts "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";
    //  input Modelica_LinearSystems2.Types method=Modelica_LinearSystems2.Types.Method.Trapezoidal
      output Modelica_LinearSystems2.DiscreteStateSpace dss(
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
      dss.Ts := Ts;
      dss.method := method;

      if method == Method.ExplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*pre(der_x)
     */
        dss.A := identity(nx) + Ts*A;
        dss.B := Ts*B;
        dss.C := C;
        dss.D := D;
        dss.B2 := zeros(nx, nu);

      elseif method == Method.ImplicitEuler then
            /*  der_x = A*x + B*u
             x = pre(x) + Ts*der_x
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          Ts*A);
        dss.B2 := LU_solve2(
              LU,
              pivots,
              Ts*B);
        dss.A := LU_solve2(
              LU,
              pivots,
              identity(nx));
        dss.B := dss.A*dss.B2;
        dss.C := C;
        dss.D := dss.C*dss.B2 + D;

      elseif method == Method.Trapezoidal then
            /*  der_x = A*x + B*u
             x = pre_x + (Ts/2)*(pre_der_x + der_x);
     */
        (LU,pivots) := Modelica_LinearSystems2.Math.Matrices.LU(identity(nx) -
          (Ts/2)*A);
        dss.B2 := LU_solve2(
              LU,
              pivots,
              (Ts/2)*B);
        dss.A := LU_solve2(
              LU,
              pivots,
              identity(nx) + (Ts/2)*A);
        dss.B := dss.A*dss.B2 + dss.B2;
        dss.C := C;
        dss.D := dss.C*dss.B2 + D;

      elseif method == Method.StepExact then
           /* x = phi*pre(x) + gamma*pre(u);
       y = C*x + D*u
    */
        (dss.A,dss.B) := Modelica.Math.Matrices.integralExp(
              A,
              B,
              Ts);
        dss.C := C;
        dss.D := D;
        dss.B2 := zeros(nx, nu);

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
        (dss.A,dss.B,dss.B2) := Modelica.Math.Matrices.integralExpT(
              A,
              B,
              Ts);
        dss.B2 := dss.B2/Ts;
        dss.B := dss.A*dss.B2 + dss.B - dss.B2;
        dss.C := C;
        dss.D := dss.C*dss.B2 + D;

      elseif method == Method.ImpulseExact then
           /* x = phi*pre(x) + phi*B*u;
        y = C*x
        (u = [1,0,1,1,0..,0])
      Limitations: The infinite impulses at t = kT is ignored in the mapping
    */

        dss.A := Modelica.Math.Matrices.exp(A, Ts);
        dss.B := dss.A*B;
        dss.C := C;
        dss.D := C*B;
        dss.B2 := B;

      else
        assert(false, "Argument method (= " + String(method) +
          ") of makeDiscrete is wrong.");
      end if;
      annotation (
         Documentation(info="<html>
<p>
This function derives a linear time invariant difference
equation system in state space form:
</p>
<blockquote><pre>
<b>xd</b>(Ts*(k+1))       = <b>Ad</b> * <b>xd</b>(Ts*k) + <b>Bd</b> * <b>ud</b>(Ts*k)
<b>yd</b>(Ts*k)           = <b>Cd</b> * <b>xd</b>(Ts*k) + <b>Dd</b> * <b>ud</b>(Ts*k)
<b>x</b>_continuous(Ts*k) =      <b>xd</b>(Ts*k) + <b>B2</b> * <b>ud</b>(Ts*k)
</pre></blockquote>
<p>
from the matrices <b>A</b>, <b>B</b>, <b>C</b>, <b>D</b> of the corresponding continuous system
</p>
<blockquote><pre>
der(<b>x</b>(t)) = <b>A</b> * <b>x</b>(t) + <b>B</b> * <b>u</b>(t)
    <b>y</b>(t)  = <b>C</b> * <b>x</b>(t) + <b>D</b> * <b>u</b>(t)
</pre></blockquote>
<p>
The function is similar to
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.'constructor'.fromStateSpace\">fromStateSpace</a>
but the inputs are restricted to the matrices, the sample time and the discretization method.
</p>

<h4>Example</h4>
<blockquote><pre>
  import dss=Modelica_LinearSystems2.DiscreteStateSpace;
  Real A[1,1] = [1];
  Real B[1,1] = [1];
  Real C[1,1] = [1];
  Real D[1,1] = [0];

public
  DiscreteStateSpace dss;

<b>algorithm</b>
  dss := dss.'constructor'.fromMatrices2(A, B, C, D);

  //or just:
  //dss := dss(A, B, C, D, Ts=0.1, method=Modelica_LinearSystems2.Types.Method.Trapezoidal);

  //  dss.A = [1.1053],
  //  dss.B = [0.11080],
  //  dss.C = [1],
  //  dss.D = [0.0526],
  //  dss.Ts = 0.1,
  //  dss.B2 = [0.0526],
  //  dss.method = Modelica_LinearSystems2.Types.Method.Trapezoidal
</pre></blockquote>
</html>"));
    end fromMatrices2;
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

  encapsulated operator '-'
    "Contains operators for subtraction of discrete state space systems"
    import Modelica;

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
<h4>Syntax</h4>
<blockquote>
<pre>
dss = DiscreteStateSpace.'-'.<b>subtract</b>(dss1, dss2)
</pre>
</blockquote>

<h4>Description</h4>
<p>
This operator function computes the subtraction of two discrete state space
systems connected in parallel, i.e. the inputs are the same and the outputs
of the two systems are subtracted. Therefore, The systems must have the same
number of inputs and outputs but not the same number of states.
The resulting system has an order of system_order1 + system_order2.
</p>
<p>
The operator is used by writing just the following command:
</p>
<blockquote><pre>
dss3 := dss1 - dss2;
</pre></blockquote>

<h4>Example</h4>
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
</html>"));
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
<p>
This package contains operators for subtraction of discrete state space records.
</p>
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
<h4>Syntax</h4>
<blockquote><pre>
    (y) = DiscreteStateSpace.<b>timeResponse</b>(dss, u)
            or
(y, xc) = DiscreteStateSpace.<b>timeResponse</b>(dss, u, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
This function calculates the time responses to input u of a discrete state space system.
Default of initial state <b>x0</b> is <b>x0</b>=<b>0</b>.
</p>

<h4>Example</h4>
<blockquote><pre>
  import dss=Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.StateSpace;
  StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
     B=[1],
     C=[2],
     D=[0]);
  Real Ts=0.1;
  dss=dss(ss,Ts);
  Real x0[1]={0};

  Real y[:,:]=dss.timeResponse(dss,ones(51,1));

//  y=[0.09524, 0.2766, 0.4408,..., 1.9844, 1.9859, 1.9872]
</pre></blockquote>
</html>",  revisions="<html>
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

  encapsulated function initialResponse
    "Compute initial response of DiscreteStateSpace system"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss "Linear system in discrete state space form";
    input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial system state";
    input Integer samples "Number of samples";
    output Real y[samples,size(dss.C, 1)]
      "System response (dimension: (input samples) x (number of outputs))";
    output Real x_continuous[samples,size(dss.A, 1)]
      "State trajectories (dimension: (input samples) x (number of states)";

  protected
    Integer i;
    Real new_x[size(dss.A, 1)];
    Real x[size(dss.A, 1)]=x0;

  algorithm
    for i in 1:samples loop
      new_x := dss.A*x;
      y[i, :] := dss.C*x;
      x_continuous[i, :] := x;
      x := new_x;
    end for;
    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
    (y) = DiscreteStateSpace.<b>initialResponse</b>(dss, x0, samples)
            or
(y, xc) = DiscreteStateSpace.<b>initialResponse</b>(dss, x0, samples)
</pre></blockquote>

<h4>Description</h4>
<p>
Function DiscreteStateSpace.initialResponse calculates the initial response to Default of initial state <b>x0</b> of a discrete state space system.
Input <b>sample</b> is the number of samples. Sample time is the sample time of the discrete state space system.
</p>

<h4>Example</h4>
<blockquote><pre>
  import dss=Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.StateSpace;
  StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1],
    B=[1],
    C=[2],
    D=[0]);
  Real Ts=0.1;
  dss=dss(ss,Ts);
  Real x0[1]={1};

  Real y[:,:]=dss.initialResponse(dss,x0,50);

//  y=[2, 1.8095, 1.6372,..., 0.01812, 0.01639, 0.01483]
</pre></blockquote>
</html>",  revisions="<html>
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

encapsulated package Analysis
    "Package of functions to analyse discrete state space system represented by a DiscreteStateSpace record"
    import Modelica;
  extends Modelica.Icons.Package;

  encapsulated function eigenValues
    "Calculate the eigenvalues of a linear discrete state space system and write them in a complex vector"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;

    input DiscreteStateSpace dss "Discrete state space system";
    output Complex eigvalues[size(dss.A, 1)] =
      Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(dss.A)
      "Eigenvalues of the system";
  algorithm

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
eigenvalues = DiscreteStateSpace.Analysis.<b>eigenValues</b>(dss)
</pre></blockquote>

<h4>Description</h4>
<p>
Calculate the eigenvalues of a discrete state space system, i.e. the eigenvalues of the system matrix <b>A</b> of a discrete state space system.
The output is a complex vector containing the eigenvalues.<br>
The eigenvalues <b>ev</b>_d of the discrete system are related to the eigenvalues of the corresponding continuous system <b>ev</b>_c by
</p>
<blockquote>
<b>ev</b>_d = exp(Ts*<b>ev</b>_c).
</blockquote>

<h4>Example</h4>
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
</html>"));
  end eigenValues;

  encapsulated function timeResponse
      "Calculate the time response of a discrete state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

    input TimeResponse response=TimeResponse.Step;
    extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete(redeclare Real
               y[
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
<h4>Syntax</h4>
<blockquote><pre>
      (y) = DiscreteStateSpace.Analysis.<b>timeResponse</b>(responseType, dss)
(y, t, x) = DiscreteStateSpace.Analysis.<b>timeResponse</b>(responseType, dss, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function timeResponse calculates the time responses of a discrete state space
system. The type of the time response is defined by the input <b>responseType</b>, i.e.
</p>
<blockquote><pre>
Impulse &quot;Impulse response&quot;,
Step    &quot;Step response&quot;,
Ramp    &quot;Ramp response&quot;,
Initial &quot;Initial condition response&quot;
</pre></blockquote>
<p>
Starting at x(t=0)=x0 and y(t=0)=C*x0 + D*u0, the outputs y and states x
are calculated for each time step t=k*dss.Ts.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A = [0.904837418036],
    B = [0.095162581964],
    C = [2],
    D = [0],
    B2 = [0],
    Ts = 0.1);

  Real tSpan = 0.4;

  Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;

  Real x0[1] = {0};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(dss,tSpan,response,x0);
//  y[:,1,1] = {0, 0.1903, 0.3625, 0.5184, 0.6594}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>
</html>"));
  end timeResponse;

  encapsulated function impulseResponse
      "Calculate the impulse time response of a discrete state space system"

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

    (y,t,x_discrete) :=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse,
          x0=zeros(size(dss.A, 1)));

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
      (y) = DiscreteStateSpace.Analysis.<b>impulseResponse</b>(dss)
(y, t, x) = DiscreteStateSpace.Analysis.<b>impulseResponse</b>(dss, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and states <b>x</b> are calculated for each time step t=k*dss.Ts. The function call
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.impulseResponse(dss, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(dss, tSpan, response=Types.TimeResponse.Impulse, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
<p>
Note that an appropriate impulse response of a discrete system that is comparable
to the impulse response of the corresponding continuous system requires
the &quot;ImpulseExact&quot; conversion from continuous system to discrete system.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A=[0.904837418036 ],
    B=[0.095162581964],
    C=[2],
    D=[0],
    B2=[0],
    Ts=0.1,
    method = Modelica_LinearSystems2.Types.Method.StepExact);

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.impulseResponse(dss,tSpan);
//  y[:,1,1]  = {0, 0.190, 0.1722, 0.1558, 0.1410}
//         t = {0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = = {0, 0.0952, 0.08611, 0.0779, 0.07050}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.impulseResponse\">StateSpace.Analysis.impulseResponse</a>
</p>
</html>"));
  end impulseResponse;

  encapsulated function stepResponse
      "Calculate the step time response of a discrete state space system"

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

    (y,t,x_discrete) :=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step,
          x0=zeros(size(dss.A, 1)));

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
      (y) = DiscreteStateSpace.Analysis.<b>stepResponse</b>(dss)
(y, t, x) = DiscreteStateSpace.Analysis.<b>stepResponse</b>(dss, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>stepResponse</b> calculates the step response of a discrete
state space system. Starting at
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and the states <b>x</b> are calculated for each
time step t=k*dss.Ts. The function call
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.stepResponse(dss, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(response=Types.TimeResponse.Step, dss, tSpan, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
<p>
Note that an appropriate step response of a discrete system that is comparable
to the step response of the corresponding continuous system requires
the &quot;StepExact&quot; conversion from continuous system to discrete system.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A=[0.904837418036 ],
    B=[0.095162581964],
    C=[2],
    D=[0],
    B2=[0],
    Ts=0.1,
    method = Modelica_LinearSystems2Types.Method.StepExact);

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x) := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.stepResponse(dss,tSpan);
//  y[:,1,1]={0, 0.19, 0.3625, 0.518, 0.659}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={0, 0.0952, 0.1813, 0.2592, 0.33}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.stepResponse\">StateSpace.Analysis.stepResponse</a>
</p>
</html>"));
  end stepResponse;

  encapsulated function rampResponse
      "Calculate the ramp time response of a discrete state space system"

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

    (y,t,x_discrete) :=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp,
          x0=zeros(size(dss.A, 1)));

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
      (y) = DiscreteStateSpace.Analysis.<b>rampResponse</b>(dss)
(y, t, x) = DiscreteStateSpace.Analysis.<b>rampResponse</b>(dss, tSpan, x0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rampResponse</b> calculates the time response
of a discrete state space system for ramp imput u = t. Starting at
<b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dss.Ts.
The function call
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.rampResponse(dss, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(response=Types.TimeResponse.Ramp, dss, tSpan, x0=fill(0,size(ss.A,1))).
</pre></blockquote>
<p>
Note that an appropriate ramp response of a discrete system that is comparable to the ramp response of the corresponding continuous system
requires the &quot;RampExact&quot; conversion from continuous system to discrete system.
</p>

<h4>Example</h4>
<blockquote><pre>
   Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
      A=[0.904837418036 ],
      B=[0.095162581964],
      C=[2],
      D=[0.0967483607192],
      B2=[0.0483741803596],
      Ts=0.1,
      method = Modelica_LinearSystems2.Types.Method.RampExact);

  Real tSpan= 0.4;

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x) := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.rampResponse(dss,tSpan);
//  y[:,1,1] = {0, 0.00967, 0.03746, 0.08164, 0.14064}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1] = {0, 0.00484, 0.01873, 0.04082, 0.07032}
</pre></blockquote>


<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.rampResponse\">StateSpace.Analysis.rampResponse</a>
</p>
</html>"));
  end rampResponse;

  encapsulated function initialResponse
      "Calculate the time response of a discrete state space system for given initial condition and zero inputs"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input Real x0[:]=fill(0, 0) "Initial state vector";

    // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask_discrete(redeclare Real
               y[
             :,size(dss.C, 1),1], redeclare Real x_discrete[:,size(dss.A, 1),1]);
    protected
    Real tSpanVar;

  algorithm
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
    else
      tSpanVar := tSpan;
    end if;

    (y,t,x_discrete) :=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial,
          x0=x0);

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
      (y) = DiscreteStateSpace.Analysis.<b>initialResponse</b>(x0, dss)
(y, t, x) = DiscreteStateSpace.Analysis.<b>initialResponse</b>(x0, dss, tSpan)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>initialResponse</b> calculates the time response of
a discrete state space system for given initial condition and zero inputs.
Starting at <b>x</b>(t=0)=<b>0</b> and <b>y</b>(t=0)=<b>C</b>*<b>x</b>0 + <b>D</b>*<b>u</b>0,
the outputs <b>y</b> and <b>x</b> are calculated for each time step t=k*dss.Ts. The function call
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.initialResponse(x0,dss, dt, tSpan)
</pre></blockquote>
<p>
gives the same result as
</p>
<blockquote><pre>
DiscreteStateSpace.Analysis.timeResponse(dss, tSpan, response=Types.TimeResponse.Initial, x0=x0).
</pre></blockquote>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A=[0.904837418036 ],
    B=[0.095162581964],
    C=[2],
    D=[0],
    B2=[0],
    Ts=0.1,
    method = Modelica_LinearSystems2.Types.Method.StepExact);

  Real tSpan= 0.4;
  Real x0[1] = {1};

  Real y[5,1,1];
  Real t[5];
  Real x[5,1,1]

<b>algorithm</b>
  (y,t,x):=Modelica_LinearSystems2.DiscreteStateSpace.Analysis.initialResponse(x0,dss,tSpan);
//  y[:,1,1]={2, 1.809, 1.637, 1.4812, 1.3402}
//         t={0, 0.1, 0.2, 0.3, 0.4}
//  x[:,1,1]={1, 0.9048, 0.8186, 0.7406, 0.6701}
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Analysis.timeResponse\">DiscreteStateSpace.Analysis.timeResponse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace.Analysis.initialResponse\">StateSpace.Analysis.initialResponse</a>
</p>
</html>"));
  end initialResponse;

end Analysis;

  encapsulated package Design
    "Package of functions to design discrete state space controllers and observers"
    import Modelica;
  extends Modelica.Icons.Package;

  encapsulated function assignPolesMI
      "Pole assignment design algorithm for multi input systems"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
  //  import Modelica.Utilities.Streams.print;
    import Modelica_LinearSystems2.TransferFunction;
    import Modelica_LinearSystems2.Math.Matrices;

    input DiscreteStateSpace dss "state space system";

    input Complex gamma[:]=fill(Complex(0), 0) "Designed Poles";
  //  input Integer np=size(gamma, 1) "number of given eigenvalues to assign";
    input Real alpha=exp(-1e10)
        "maximum admissible value for the moduli of the eigenvalues of A which will not be modified by the eigenvalue assignment algorithm";
    input Real tolerance=Modelica.Math.Matrices.norm(dss.A, 1)*1e-12
        "The tolerance to be used in determining the controllability of (A,B)";
    input Boolean calculateEigenvectors=false
        "Calculate the eigenvectors X of the closed loop system when true";

    output Real K[size(dss.B, 2),size(dss.A, 1)]
        "State feedback matrix assigning the desired poles";
    output Real S[:,:] "Closed loop System matrix";
    output Complex po[size(dss.A, 1)] "poles of the closed loop system";
    output Integer nfp
        "number of eigenvalues that are not modified with respect to alpha";
    output Integer nap "number of assigned eigenvalues";
    output Integer nup "number of uncontrollable eigenvalues";
    output Complex X[size(dss.A, 1),size(dss.A, 1)]
        "eigenvectors of the closed loop system";

    protected
    Real A_rsf[size(dss.A, 1),size(dss.A, 2)];
    Real B_rsf[size(dss.B, 1),size(dss.B, 2)];
    Real Q[size(dss.A, 1),size(dss.A, 1)];
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
    Integer n=size(dss.A, 1);
    Integer nccA "number of conjugated complex pole pairs of openloop system";
    Integer nccg "number of conjugated complex pole pairs of gamma";
    Integer rpg "number of real poles in gamma";
    Integer rpA "number of real poles of open loop system";
    Integer ncc "Min(nccA, nccg)";
    Integer rp "Min(rpg, rpA)";
    Integer ng=size(gamma,1);
    Integer nr "Differenz between rpA and rpg; Sign(rpA-rpg)*(rpA-rpg)";

    Real alphaReal[size(dss.A, 1)]
        "Real part of eigenvalue=alphaReal+i*alphaImag";
    Real alphaImag[size(dss.A, 1)]
        "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

    Complex SS[:,:];
    Complex Xj[:,:];
    Complex h;

    Real dist;
    Real evImag;

  algorithm
    assert(size(gamma, 1) <= size(dss.A, 1),
      "At most n (order of dss) eigenvalues can be assigned");

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

    // put matrix dss.A to real Schur form A <- QAQ' and compute B <- QB
    (A_rsf,Z,alphaReal,alphaImag) := Matrices.rsf2(dss.A);
    ZT := transpose(Z);

    // reorder real Schur form according to alpha
    (A_rsf,Z,alphaReal,alphaImag) := Matrices.Internal.reorderRSFd(
        A_rsf,
        identity(size(A_rsf, 1)),
        alphaReal,
        alphaImag,
        alpha);
    ZT := transpose(Z)*ZT;
    B_rsf := ZT*dss.B;

    // determine number of poles not to be assigned according to alpha
    nfp := 0;
    for i in 1:n loop
      if alphaReal[i]^2+alphaImag[i]^2 < alpha^2 then
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
          identity(n - nfp));//The Schur vector matrix is identity, since A_rsf already has Schur form

      A_rsf[1:nfp, nfp + 1:n] := A_rsf[1:nfp, nfp + 1:n]*Q2;
      B_rsf[nfp + 1:n, :] := transpose(Q2)*B_rsf[nfp + 1:n, :];
      ZT[nfp + 1:n, :] := transpose(Q2)*ZT[nfp + 1:n, :];
    end if;

    // main algorithm
    K := zeros(size(dss.B, 2), size(dss.A, 1));
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

        Ks1 := DiscreteStateSpace.Internal.assignOneOrTwoPoles(
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

        Ks2 := DiscreteStateSpace.Internal.assignOneOrTwoPoles(
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

        Ks2 := DiscreteStateSpace.Internal.assignOneOrTwoPoles(
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
  //      Modelica.Utilities.Streams.print("counter2Case3 = " + String(counter2));
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

      Ks2 := DiscreteStateSpace.Internal.assignOneOrTwoPoles(
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

    S := dss.A - dss.B*K;
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
(K, S, po, nfp, nap, nup) = DiscreteStateSpace.Design.<b>assignPolesMI</b>(dss, gamma, np, tol, calculateEigenvectors)
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
At the beginning of the algorithm, the feedback matrix <b>K</b> is set
to zero (<b>K</b> = <b>0</b>) and the matrix <b>A</b> is reduced to an ordered
real Schur form by separating its spectrum in two parts
</p>
<blockquote><pre>
              | <b>F</b>1  <b>F</b>3|
 <b>F</b> = <b>Q</b>*<b>A</b>*<b>Q</b>' = |       |
              | <b>0</b>   <b>F</b>2|
</pre>
</blockquote>
<p>
in such a way, that <b>F</b>1 contains the eigenvalues that will be retained
and <b>F</b>3 contains the eigenvalues going to be modified. On the suggestion
of [1] the eigenvalues <i>evr</i> to be retained are chosen as
</p>
<blockquote><pre>
evr = {s in C: Re(s) &lt; -alpha, alpha &gt;=0}
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
appropriate feedback matrix <b>K</b>i. The matrix <b>F</b> - <b>G</b>*<b>K</b>i remains
in real Schur form. The assigned eigenvalue(s) is (are) then moved to another diagonal
position of the real Schur form using reordering techniques <b>F</b>
&lt; -- <b>Q</b>i*<b>F</b>*<b>Q</b>i'  and a new block is transferred to the
lower right diagonal position. The transformations are accumulated in <b>Q</b>i
and are also applicated to the matrices
</p>
<blockquote><pre>
<b>G</b> &lt; - <b>Q</b>i*<b>G</b> <b>Q</b> &lt; - <b>Q</b>i*<b>Q</b>
</pre></blockquote>
<p>
The eigenvalue(s) to be assigned at  each step is (are) chosen such that
the norm of each <b>K</b>i is minimized [1].
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A=[-1, 1, 1;0, 1, 1;0, 0, 1],
    B=[0; 0; 1],
    C=[0, 1, 0],
    D=[0]);

  Real Q[3,3];

<b>algorithm</b>
  Q := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.observabilityMatrix(dss);
// Q = [0, 1, 0; 0, 1, 1; 1, 1, 2]
</pre></blockquote>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga A. (1981):</dt>
<dd> <b>A Schur method for pole assignment</b>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>
</html>",revisions="<html>
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

  end Design;

encapsulated package Plot
    "Package of functions to plot discrete state space system responses"
    import Modelica;
  extends Modelica.Icons.Package;
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
        "True, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min=0.1
        "Minimum frequency value, if autoRange = false";
  input Modelica.SIunits.Frequency f_max=10
        "Maximum frequency value, if autoRange = false";

  input Boolean magnitude=true "= true, to plot the magnitude of dtf"
                                                                     annotation(choices(checkBox=true));
  input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot());

  input Boolean Hz=true
        "= true, to plot abszissa in [Hz], otherwise in [rad/s] (= 2*pi*Hz)"
                                                                         annotation(choices(checkBox=true));
  input Boolean dB=false
        "= true, to plot magnitude in [], otherwise in [dB] (=20*log10(value))"
                                                                            annotation(choices(checkBox=true),Dialog(enable=magnitude));
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
    Hz=Hz,
    magnitude=magnitude,
    dB=dB,
    phase=phase,
    defaultDiagram=defaultDiagram,
    device=device);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>bodeSISO</b>(dss)
   or
DiscreteStateSpace.Plot.<b>bodeSISO</b>(
  dss,
  iu,
  iy,
  nPoints,
  autoRange,
  f_min,
  f_max,
  magnitude=true,
  phase=true,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a bode-diagram of the transfer function corresponding
to the behavior of the state space system from iu'th element of the input
vector <b>u</b> to the iy'th element of the output vector <b>y</b>.
</p>

<h4>Example</h4>
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
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodeMagDis.png\">
</p>
<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/bodePhaseDis.png\">
</p>
</html>"));
end bodeSISO;

encapsulated function timeResponse
      "Plot the time response of a discrete state space system. The response type is selectable"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteStateSpace dss;
  input Real tSpan=0 "Simulation time span [s]";

  input TimeResponse response=TimeResponse.Step;

  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

  input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(choices(checkBox=true));

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

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>timeResponse</b>(dss);
   or
DiscreteStateSpace.Plot.<b>timeResponse</b>(
  dss,
  tSpan,
  response,
  x0,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the time response of a discrete state space system. The character of the time response if defined by the input <tt>response</tt>, i.e. Impulse, Step, Ramp, or Initial.
</p>

<h4>Example</h4>
<blockquote><pre>
  import Modelica_LinearSystems2.DiscreteStateSpace;

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
  DiscreteStateSpace.Plot.timeResponse(dss, response=response);
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
end timeResponse;

  encapsulated function impulse
      "Impulse response plot of a discrete state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";
    input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

    input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Impulse response"));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Impulse "type of time response";
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

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>impulse</b>(dss);
   or
DiscreteStateSpace.Plot.<b>impulse</b>(
  dss,
  tSpan,
  x0,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the impulse responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
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
  Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse(dss)
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end impulse;

  encapsulated function step
      "Step response plot of a discrete state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";

    input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Step response"));

    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "type of time response";

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

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>step</b>(dss);
   or
DiscreteStateSpace.Plot.<b>step</b>(
  dss,
  tSpan,
  x0,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the discrete step responses of a state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
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
  Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss, tSpan=3)
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initialResponse\">initialResponse</a>
</p>
</html>"));
  end step;

  encapsulated function ramp
      "Ramp response plot of a discrete state space system"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;
      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteStateSpace dss;
    input Real tSpan=0 "Simulation time span [s]";

    input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(choices(checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
           heading="Ramp response"));
    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Ramp "type of time response";

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

    annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>ramp</b>(ss);
   or
DiscreteStateSpace.Plot.<b>ramp</b>(
  dss,
  tSpan,
  x0,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the ramp responses of a discrete state space system for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>. See also
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a> and
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.initialResponse\">initialResponse</a>.
</p>

<h4>Example</h4>
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
</html>"));
  end ramp;

encapsulated function initialResponse
      "Initial condition response plot of a discrete state space system"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.Utilities.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

  input DiscreteStateSpace dss;
  input Real tSpan=0 "Simulation time span [s]";
  input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1)) "Initial state vector";

  input Boolean subPlots=true
        "True, if all subsystem time responses are plotted in one window with subplots"
                                                                                     annotation(choices(checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunctionMIMO(defaultDiagram=Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(
        heading="Initial response"));
    protected
      Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Initial "type of time response";

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

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
DiscreteStateSpace.Plot.<b>initialResponse</b>(ss);
   or
DiscreteStateSpace.Plot.<b>initialResponse</b>(
  dss,
  tSpan,
  x0,
  defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramPolesAndZeros\">Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse</a>(),
  device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>())
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots the initial responses of a discrete state space system for the initial state vector x0 for each system corresponding to the transition matrix. It is based on <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.timeResponse\">timeResponse</a>.
</p>

<h4>Example</h4>
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

<h4>See also</h4>
<p>
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse\">impulse</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.step\">step</a>,
<a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp\">ramp</a>
</p>
</html>"));
end initialResponse;

end Plot;

  encapsulated package Conversion
    "Package of functions for conversion of DiscreteStateSpace data record"
    import Modelica;
  extends Modelica.Icons.Package;

    encapsulated function toDiscreteZerosAndPoles
      "Generate a discrete zeros-and-poles representation from a discrete SISO state space representation"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.ZerosAndPoles;
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

        poles :=Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.eigenValues_dhseqr(ssm.A);
          //ssm.A is of upper Hessenberg form
        zeros := StateSpace.Internal.invariantZeros2(ssm);
        cpoles := fill(Complex(0),size(poles,1));
        czeros := fill(Complex(0),size(zeros,1));

        if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
          assert(size(ss.B, 2) == 1, "Function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
             + String(size(ss.B, 2)) + " instead of 1");
          assert(size(ss.C, 1) == 1, "Function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
             + String(size(ss.C, 1)) + " instead of 1");
        end if;
        dzp := DiscreteZerosAndPoles(
            z=zeros,
            p=poles,
            k=1,
            Ts=dss.Ts, method=dss.method);
    // set frequency to a complex value which is whether pole nor zero
        for i in 1:size(poles,1) loop
          cpoles[i] := if Modelica.ComplexMath.'abs'(poles[i]) > 0 then Modelica.ComplexMath.log(poles[i])/dss.Ts else Complex(-100);
        end for;
        for i in 1:size(zeros,1) loop
          czeros[i] := if Modelica.ComplexMath.'abs'(zeros[i]) > 0 then Modelica.ComplexMath.log(zeros[i])/dss.Ts else Complex(-100);
        end for;

         v := sum(cat(1, czeros[:].re,  cpoles[:].re))/max(size(czeros,1)+size(cpoles,1),1) + 13/19;
    //     v := sum(cat(1, zeros[:].re,  poles[:].re))/max(size(zeros,1)+size(poles,1),1);
        frequency := Complex(v)*17/19;
        cfrequency :=Modelica.ComplexMath.exp(frequency*dss.Ts);
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

      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteStateSpace.Conversion.<b>toDiscreteZerosAndPoles</b>(dss)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a DiscreteZerosAndPoles record
</p>
<blockquote><pre>
           product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
dzp = k * ---------------------------------------------------------
           product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre></blockquote>
<p>
of a system from discrete state space representation using the transformation algorithm described in [1].
</p>
<p>
The uncontrollable and unobservable parts are isolated and the eigenvalues and invariant zeros of the controllable and observable sub system are calculated.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(
    A = [0.9048, 0.0,    0.0;
         0.0,    0.8187, 0.0;
         0.0,    0.0,    0.7408],
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

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>
</html>",     revisions="<html>
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
    end toDiscreteZerosAndPoles;

    encapsulated function toDiscreteZerosAndPolesMIMO
      "Generate a zeros-and-poles representation from a MIMO state space representation"

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;
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
      annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
dzp = DiscreteStateSpace.Conversion.<b>toDiscreteZerosAndPolesMIMO</b>(dss)
</pre></blockquote>

<h4>Description</h4>
<p>
Computes a matrix of DiscreteZerosAndPoles records
</p>
<blockquote><pre>
           product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
dzp = k * ---------------------------------------------------------
           product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre></blockquote>
<p>
of a system from discrete state space representation, i.e. isolating the uncontrollable and unobservable parts and the eigenvalues and invariant zeros of the controllable and observable sub systems are calculated. The algorithm applies the method described in [1] for each single-input-output pair.
</p>

<h4>Example</h4>
<blockquote><pre>
  Modelica_LinearSystems2.DiscreteStateSpace ss=Modelica_LinearSystems2.DiscreteStateSpace(
    A = [0.9048, 0.0,    0.0;
         0.0,    0.8187, 0.0;
         0.0,    0.0,    0.7408]

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

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
</dl>
</html>",     revisions="<html>
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
    end toDiscreteZerosAndPolesMIMO;

  function toDiscreteTransferFunction
      "Generate a TransferFunction data record from a SISO DiscreteStateSpace data record"

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
<h4>Syntax</h4>
<table>
<tr> <td align=right>  dtf </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteTransferFunction</b>(dss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a DiscreteTransferFunction record
<blockquote><pre>
           n(z)     b0 + b1*z + ... + bn*z^n
  dtf = -------- = --------------------------
           d(z)     a0 + a1*z + ... + an*z^n
 </pre></blockquote>

The algorithm uses <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles\">toDiscreteZerosAndPoles</a> to convert the
discrete state space system into a discrete zeros and poles representation first and after that <a href=\"Modelica://Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction\">DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction</a> to generate the transfer function.



<h4>Example</h4>
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
<h4>Syntax</h4>
<table>
<tr> <td align=right>  dtf </td><td align=center> =  </td>  <td> DiscreteStateSpace.Conversion.<b>toDiscreteTransferFunctionMIMO</b>(dss)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a matrix of DiscreteTransferFunction records
<blockquote><pre>
           n(z)     b0 + b1*z + ... + bn*z^n
  dtf = -------- = --------------------------
           d(z)     a0 + a1*z + ... + an*z^n
 </pre></blockquote>
with repetitive application of <a href=\"Modelica://Modelica_LinearSystems2.DiscreteStateSpace.Conversion.toDiscreteTransferFunction\">Conversion.toDiscreteTransferFunction</a>


<h4>Example</h4>
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
  end toDiscreteTransferFunctionMIMO;

  end Conversion;

  encapsulated package Import
    "Package of functions to generate a DiscreteStateSpace data record from imported data"
    import Modelica;
    extends Modelica.Icons.Package;

  encapsulated function fromFile
    "Read a DiscreteStateSpace data record from mat-file"
    import Modelica;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteStateSpace;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension;

    input String fileName = "dslin.mat" "Name of the state space system data file"
      annotation (
        Dialog(
          loadSelector(
            filter="MAT files (*.mat);; All files (*.*)",
            caption="State space system data file")));
    input String matrixName = "ABCD"
      "Name of the state space system matrix" annotation(Dialog);
    protected
    Integer xuy[3] = ReadSystemDimension(fileName, matrixName) annotation(__Dymola_allowForSize=true);
    Integer nx = xuy[1] annotation(__Dymola_allowForSize=true);
    Integer nu = xuy[2] annotation(__Dymola_allowForSize=true);
    Integer ny = xuy[3] annotation(__Dymola_allowForSize=true);

    public
    output DiscreteStateSpace result(
      redeclare Real A[nx, nx],
      redeclare Real B[nx, nu],
      redeclare Real B2[nx, nu],
      redeclare Real C[ny, nx],
      redeclare Real D[ny, nu]) "Model read from file";

    protected
    Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(
      fileName, matrixName, nx + ny, nx + nu);
    Real B2[nx,nu] = Streams.readRealMatrix(fileName, "B2", nx, nu);
    Real Ts[1,1] = Streams.readRealMatrix(fileName, "Ts", 1, 1);

  algorithm
    result.A := ABCD[1:nx, 1:nx];
    result.B := ABCD[1:nx, nx + 1:nx + nu];
    result.B := B2;
    result.C := ABCD[nx + 1:nx + ny, 1:nx];
    result.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    result.Ts := scalar(Ts);
    Streams.print("StateSpace record loaded from file: \"" +
      Modelica.Utilities.Files.fullPathName(fileName) + "\"");

    annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  dss </td><td align=center> =  </td>  <td> DiscreteStateSpace.Import.<b>fromFile</b>(fileName, matrixName)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Reads and loads a discrete state space system from a mat-file <tt>fileName</tt>. The discrete System on the file must be described as:
</p>

<blockquote><pre>
<b>x</b>(Ts*(k+1))  = <b>A</b> * <b>x</b>(Ts*k) + <b>B</b> * <b>u</b>(Ts*k)
<b>y</b>(Ts*k)         = <b>C</b> * <b>x</b>(Ts*k) + <b>D</b> * <b>u</b>(Ts*k)
<b>x</b>_continuous(Ts*k) =  <b>x</b>(Ts*k) + <b>B2</b> * <b>u</b>(Ts*k)
</pre></blockquote>
<p>
with (for more details see <a href=\"modelica://Modelica_LinearSystems2.DiscreteStateSpace.'constructor'.fromStateSpace\">DiscreteStateSpace.'constructor'.fromStateSpace</a>):
</p>
<ul>
<li> <b>Ts</b> - the sample time,</li>
<li> <b>k</b> - the index of the actual sample instance (k=0,1,2,3,...),</li>
<li> <b>t</b> - the time,</li>
<li> <b>u</b>(t) - the input vector,</li>
<li> <b>y</b>(t) - the output vector,</li>
<li> <b>x</b>(t) - the discrete state vector (x(t=Ts*0) is the initial state),</li>
<li> <b>x</b>_continuous(t) - the state vector of the continuous system
     from which the discrete block has been derived,</li>
<li> <b>A, B, C, D, B2</b> - matrices of appropriate dimensions.</li>
</ul>

<p>
The file must contain
</p>
<ul>
<li> the Real matrix [A, B; C, D]  with name &quot;matrixName&quot;,</li>
<li> the Integer matrix &quot;nx[1,1]&quot; defining the number of states (that is the number of rows of the square matrix A),</li>
<li> the Real matrix B2 that has the same dimensions as B,</li>
<li> the Real matrix &quot;Ts[1,1]&quot; defining the sample time in [s] with which the continuous-time system was discretized to arrive at this discrete system</li>
</ul>

<h4>Example</h4>
<blockquote><pre>
<b>algorithm</b>
  file := Modelica.Utilities.Files.loadResource(
    &quot;modelica://Modelica_LinearSystems2/Resources/Data/dss.mat&quot;)
  dss := Modelica_LinearSystems2.DiscreteStateSpace.Import.fromFile(file)
//  dss=StateSpace(
      A=[-4.5, 1.5, 4.0; -4.0, 1.0, 4.0; -1.5, -0.5, 1],
      B=[2; 1; 2],
      C=[1, 0, 0],
      D=[0],
      B2=[0;0;0],
      Ts=0.2,
      method=Trapezoidal);
</pre></blockquote>
</html>"));
  end fromFile;

    function fromModel
      "Generate a DiscreteStateSpace data record by linearization of a modelica model"

      import Modelica;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input String modelName "Name of the Modelica model"
        annotation (Dialog(__Dymola_translatedModel(translate=true)));
      input Modelica.SIunits.Time T_linearize=0
        "point in time of simulation to linearize the model";
      input String fileName="dslin" "Name of the result file";
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.Trapezoidal "Discretization method";

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

      Integer xuy[3] = Modelica_LinearSystems2.Internal.Streams.ReadSystemDimension(fileName2, "ABCD");
      Integer nx = xuy[1];
      Integer nu = xuy[2];
      Integer ny = xuy[3];
      Real ABCD[nx + ny,nx + nu] = Streams.readRealMatrix(
        fileName2, "ABCD", nx + ny, nx + nu);
      String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);

      StateSpace ss(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,nu],
        redeclare Real C[ny,nx],
        redeclare Real D[ny,nu]) "Model linearized at initial point";
    public
      output DiscreteStateSpace result(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,nu],
        redeclare Real B2[nx,nu],
        redeclare Real C[ny,nx],
        redeclare Real D[ny,nu]) "Discrete model linearized at initial point";

    algorithm
      ss.A := ABCD[1:nx, 1:nx];
      ss.B := ABCD[1:nx, nx + 1:nx + nu];
      ss.C := ABCD[nx + 1:nx + ny, 1:nx];
      ss.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
      ss.uNames := xuyName[nx + 1:nx + nu];
      ss.yNames := xuyName[nx + nu + 1:nx + nu + ny];
      ss.xNames := xuyName[1:nx];

      result := DiscreteStateSpace(
        ss=ss,
        Ts=Ts,
        method=method);

      annotation (__Dymola_interactive=true, Documentation(info=
                                                   "<html>
<h4>Syntax</h4>
<blockquote><pre>
dss = DiscreteStateSpace.Import.<b>fromModel</b>(modelName, T_linearize, fileName, Ts, method)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate a discrete state space data record by linearization of
a model defined by modelName. The linearization is performed at time T_linearize
of the simulation. The result of linearization is transformed into a state space
record and then converted into a discrete state space record.
</p>

<h4>Example</h4>
<blockquote><pre>
  String modelName = &quot;Modelica_LinearSystems2.Utilities.Plants.DoublePendulum&quot;;
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
</html>"));
    end fromModel;

  end Import;

  encapsulated package Internal
    "Package of internal material of record DiscreteStateSpace (for advanced users only)"
    extends Modelica.Icons.InternalPackage;
    import Modelica;

    function timeResponseSamples
      "Estimate reasonable discretisation sample time and simulation time span for time response plot"
      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;

      input Modelica_LinearSystems2.DiscreteStateSpace dss;
      output Real tSpan "Time span";
    protected
      Complex eig[size(dss.A, 1)];
      Real realp[size(dss.A, 1)];
      Real sorted[size(dss.A, 1)];
      Real indices[size(dss.A, 1)];
      Integer i;
    algorithm
      eig := Modelica_LinearSystems2.DiscreteStateSpace.Analysis.eigenValues(dss);
      for i in 1:size(dss.A, 1) loop
        eig[i] :=if Modelica.ComplexMath.'abs'(eig[i]) > 1e-10 then Modelica.ComplexMath.log(eig[i])/dss.Ts else Complex(-100);
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
      input Real x0[size(dss.A, 1)]=zeros(size(dss.A, 1))
        "Initial system state";
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

  encapsulated function assignOneOrTwoPoles
    "Algorithm to assign p (p = 1 or 2) eigenvalues"

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2;
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
      "\n In function DiscreteStateSpace.Internal.assignOneOrTwoPoles() matrix F is of size ["
       + String(size(F, 1)) + "," + String(size(F, 1)) + "] and " + String(
      size(F, 1)) + " demanded assigned poles are expected. However, " +
      String(size(gamma, 1)) + " poles are given");
  //assert(not Modelica.Math.Matrices.isEqual(G,zeros(size(G,1),size(G,2)),tolerance),"A subsystem (F, G) in DiscreteStateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix ");
    if size(gamma, 1) == 1 then
      assert(gamma[1].im == 0, "\n In function DiscreteStateSpace.Internal.assignOneOrTwoPoles() matrix F has size [" + String(size(F, 1)) + "," + String(size(F, 1)) +
        "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
         + String(gamma[1].im));
    elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
      assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im,
        "\nThe assigned pole pair given in function DiscreteStateSpace.Internal.assignOneOrTwoPoles() must be conjungated complex. However, the poles are\npole1 = "
         + String(gamma[1]) + "\npole2 = " + String(gamma[2]) +
        ". \nTry\npole1 = " + String(gamma[1]) + "\npole2 = " + String(
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
        tolerance), "A subsystem in DiscreteStateSpace.Internal.assignOneOrTwoPoles() is not controllable");

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
        Modelica.Utilities.Streams.print("\n A subsystem (F, G) in DiscreteStateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, "   + String(F[1, 1]) + " remains and " + String(gamma[1].re) + " cannot be realized");
      else
        system_ev :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(F);
        Modelica.Utilities.Streams.print("\n A subsystem (F, G) in DiscreteStateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
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
  end Internal;

  annotation (
    Icon(
      graphics={
        Rectangle(
          lineColor={160,160,164},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid,
          extent={{-100,-100},{100,100}},
          radius=25),
        Text(
          lineColor={255,255,170},
          extent={{-90,-50},{90,50}},
          textString="ss")}));
end DiscreteStateSpace;
