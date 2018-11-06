within Modelica_LinearSystems2.Examples.StateSpace;
function plotZeros "Case studies of systems with zeros"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.StateSpace;
  import Complex;

protected
  parameter Real sampleT=0.001;
  TransferFunction tf=TransferFunction({1,-0.1,1.0025}, {1,3,3,1});
  StateSpace ss=StateSpace(tf);
  Complex invZeros[:]=Modelica_LinearSystems2.StateSpace.Analysis.invariantZeros( ss);
  Complex invZero1;
  Real invMat[2*size(ss.A, 1),2*size(ss.A, 1)];
  Real t[:]=0:sampleT:20 "Time vector: (number of samples)";
  Real u[size(t, 1),1];
  Real u2[size(t, 1),1];
  Real y[:,1];
  Real y2[:,1];
  Integer windowID[size(ss.C, 1)] "ID of the used plot window";
  Real x0[size(ss.A, 1)]=zeros(size(ss.A, 1)) "Initial state vector";
  Real x[size(t,1),size(ss.A,1)];
  Real Ts=sampleT;
  Modelica_LinearSystems2.DiscreteStateSpace sd=Modelica_LinearSystems2.DiscreteStateSpace(ss, Ts, method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact);
  Boolean ok;

algorithm
  Modelica.Utilities.Streams.print(String(ss));
  Modelica.Utilities.Streams.print(String(tf));
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("Invariant zeros", invZeros);

  invZero1 := invZeros[1];
  invMat := Modelica.Math.Matrices.inv([invZero1.re*identity(size(ss.A, 1)) -
    ss.A,invZero1.im*identity(size(ss.A, 1)); -invZero1.im*identity(size(ss.A,
    1)),invZero1.re*identity(size(ss.A, 1)) - ss.A]);
  x0 := vector(2*(invMat[1:size(ss.A, 1), 1:size(ss.A, 1)])*ss.B);
  Modelica_LinearSystems2.Math.Vectors.printVector(
    x0,
    6,
    "x0");

  u[:, 1] := 2*exp(invZero1.re*t).*vector(cos(t));
  (y,x) := Modelica_LinearSystems2.DiscreteStateSpace.timeResponse(
    sd,
    u,
    x0);

    Modelica_LinearSystems2.Utilities.Plot.diagram(
        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
                  curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                           x=t,
                           y=y[:,1],
                           legend="y"),
                         Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                           x=t,
                           y=u[:,1],
                           legend="u")},
                  heading="y and u",
                  xLabel="t [s]",
                  yLabel="y and u"));

 Modelica_LinearSystems2.Utilities.Plot.diagram(
       Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
                 curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=x[:,1],
                          legend="x1"),
                        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=x[:,2],
                          legend="x2"),
                        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=x[:,3],
                          legend="x3")},
                 heading="x",
                 xLabel="t [s]",
                 yLabel="x"));

  u[:, 1] := 2*exp(invZero1.re*t).*vector(cos(2*t));

  y := Modelica_LinearSystems2.DiscreteStateSpace.timeResponse(
    sd,
    u,
    x0);

  Modelica_LinearSystems2.Utilities.Plot.diagram(
       Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
                 curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=y[:,1],
                          legend="y"),
                        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=u[:,1],
                          legend="u")},
                 heading="y and u",
                 xLabel="t [s]",
                 yLabel="y, u"));

 u[:, 1] := 2*exp(invZero1.re*t).*vector(cos(0.5*t));

  y := Modelica_LinearSystems2.DiscreteStateSpace.timeResponse(
    sd,
    u,
    x0);

   Modelica_LinearSystems2.Utilities.Plot.diagram(
       Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
                 curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=y[:,1],
                          legend="y"),
                        Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                          x=t,
                          y=u[:,1],
                          legend="u")},
                 heading="y and u",
                 xLabel="t [s]",
                 yLabel="y, u"));

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
Computes the initial condition response of the system
StateSpace <i>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</i> to the initial condition <i>x0=[1;1]</i>.
</p>

<p>
This example plots the output y and the states (x1, x2, x3) of a system with the input
</p>
<blockquote><pre>
u(t) = uk*exp(zk*t)
</pre></blockquote>
<p>
where zk is an invariant zero of the system. Assuming appropriate initial conditions, the output of the system is forced to zero. It is demonstrated that the output can also be forced to zero by applying a transient unstable input. Although the output is zero, the states show transient and unstable behavior. In comparison, the outputs as an reaction of inputs with half or double frequency are not equal to zero.
</p>
</html>"));
end plotZeros;
