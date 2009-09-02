within Modelica_LinearSystems2.Internal;
record StateSpaceR
  "Continuous state space description of a linear, time invariant differential equation system (data + operations)"

  extends Modelica.Icons.Record;

  Real A[:,:]   annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real B[size(A, 1),:]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real C[:,size(A, 1)]  annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Real D[size(C, 1),size(B, 2)] annotation(Dialog(group="der(x) = A*x + B*u;  y = C*x + D*u"));
  Integer r=0;

  annotation (defaultComponentName="stateSpaceR", Documentation(info="<html>
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
Compared to record StateSpace StateSpaceR cantains an additional Integer value r
</html>"));

  encapsulated function constructor
    "Default constructor for a StateSpace record"
    import Modelica;
    import Modelica_LinearSystems2.Internal.StateSpaceR;

    input Real A[:,size(A, 1)] "der(x) = A*x + B*u";
    input Real B[size(A, 1),:] "der(x) = A*x + B*u";
    input Real C[:,size(A, 1)] "y = C*x + D*u";
    input Real D[size(C, 1),size(B, 2)] "y = C*x + D*u";
    input Integer r;
    output StateSpaceR result(
      redeclare Real A[size(A, 1),size(A, 2)],
      redeclare Real B[size(B, 1),size(B, 2)],
      redeclare Real C[size(C, 1),size(C, 2)],
      redeclare Real D[size(D, 1),size(D, 2)]);
  algorithm
    result.A := A;
    result.B := B;
    result.C := C;
    result.D := D;
    result.r := r;
  end constructor;
end StateSpaceR;
