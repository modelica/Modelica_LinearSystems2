within Modelica_LinearSystems2.Controller;
block ZerosAndPoles
  "Continuous or discretized, single input single output block described by a ZerosAndPoles object"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Utilities.Types;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Controller.Internal;
  import Complex;

  extends Modelica_LinearSystems2.Controller.Interfaces.PartialSISO2(
      discretePart(
      x_start=x_start,
      y_start={y_start},
      ABCD=ZerosAndPoles.Conversion.toMatrices(system)));
  parameter ZerosAndPoles system "Data defining the ZerosAndPoles object";
  final parameter Integer nx=size(system.d1,1) + 2*size(system.d2,1)
    "Number of states";
  parameter Real x_start[nx] = zeros(nx) "Initial or guess values of states"
    annotation(Dialog(tab="Advanced options"));
  parameter Real y_start = 0 "Initial or guess values of output"
    annotation(Dialog(tab="Advanced options"));
  Modelica.Blocks.Interfaces.RealOutput x[nx] "State of continuous block";

protected
parameter Boolean withDelay=false;
  parameter Integer n_num1 = size(system.n1,1);
  parameter Integer n_num2 = size(system.n2,1);
  parameter Integer n_den1 = size(system.d1,1);
  parameter Integer n_den2 = size(system.d2,1);
  parameter Integer n_num = n_num1 + 2*n_num2;
  parameter Integer n_den = n_den1 + 2*n_den2;
  parameter Integer i_d = if n_num2 > n_den2 then 2*(n_num2 - n_den2) + 1 else 1;
  parameter Integer i_k = if n_num2 > n_den2 then n_den2 - (n_num2-n_den2) else n_den2;
  parameter Real eps=1e-6;
  parameter Real num[nx,2] = [system.n2;
                              [system.n1, zeros(n_num1)];
                              zeros(nx-n_num2-n_num1,2)]
    "Numerator matrix, in order that indices are defined in all situations in all if clauses";
  parameter Real den[nx,2] = [system.d2;
                              [system.d1, zeros(n_den1)];
                              zeros(nx-n_den2-n_den1,2)]
    "Denominator matrix, in order that indices are defined in all situations in all if clauses";
  Real uu[i_k + n_den1 + 1]
    "Input signals into the connected first and second order blocks";
  parameter Real k[size(uu,1)-1](each fixed = false)
    "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
  parameter Real k_total=system.k/product(k);

initial equation
  /* Compute scaling factor for every block in order that the gain of the block is 1.
     The for blocks and the if-blocks have the same structure as in the
     equation part below
  */
   for i in 1:max(n_den2, n_num2) loop
      // State space systems of order 2
      if i <= n_den2 then
        if i <= n_num2 then
            // State space system in form (1)
          k[i] =  Internal.scaleFactor2(
              num[i, 1],
              num[i, 2],
              den[i, 1],
              den[i, 2],eps);
        elseif 2*(i - n_num2) <= n_num1 then
            // State space system in form (1) with 2 first order numerator polynomials
          k[i] =  Internal.scaleFactor2(
              num[max(1,2*(i - n_num2)-1), 1] + num[max(1,min(nx,2*(i - n_num2))), 1],
              num[max(1,2*(i - n_num2)-1), 1]*num[max(1,min(nx,2*(i - n_num2))), 1],
              den[i, 1],
              den[i, 2],eps);
        elseif  2*(i-n_num2) -1== n_num1 then
            // State space system in form (2) with 1 first order numerator polynomial
          k[i] =  Internal.scaleFactor2(
              1,
              num[max(1,min(nx,2*i-n_num2-1)), 1],
              den[i, 1],
              den[i, 2],eps);
        else
            // State space system in form (3)
          k[i] =  Internal.scaleFactor2(
              1,
              1,
              den[i, 1],
              den[i, 2],eps);
        end if;
      else
         // State space system in form (1) with 2 first order denominator polynomials
        k[i] =  Internal.scaleFactor2(
            num[i, 1],
            num[i, 2],
            den[max(1,2*(i - n_den2)-1), 1] + den[max(1,2*(i - n_den2)), 1],
            den[max(1,2*(i - n_den2)-1), 1]*den[max(1,2*(i - n_den2)), 1],eps);
      end if;
    end for;

    for i in i_d:n_den1 loop
      // State space systems of order 1
      if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
         // State space system in form (4)
        k[i_k + i] =  Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1], den[n_den2 + i, 1],eps);
      elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
         // State space system in form (4)
        k[i_k + i] =  Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
          1], den[n_den2 + i, 1],eps);
      else
         // State space system in form (5)
        k[i_k + i] =  Internal.scaleFactor1(1, den[n_den2 + i, 1],eps);
      end if;
    end for;

equation
  assert(n_num <= n_den, "ZerosAndPoles transfer function is not proper as required from StateSpace system:\n"
                         + "  numerator degree (= " + String(n_num) + ") <= denominator degree (= "
                         + String(n_den) +") required.");
  if continuous then
     for i in 1:max(n_den2,n_num2) loop
        // Construct state space systems of order 2
        der(x[2*i-1]) = x[2*i];
        if i <= n_den2 then

           if i <= n_num2 then
              // State space system in form (1)
             der(x[2*i]) = if abs(den[i, 2])>eps and abs(num[i, 2])>eps then den[i, 2]*uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i] else uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i];
              uu[i+1] = if abs(den[i, 2])>eps and abs(num[i, 2])>eps then k[i]*(((num[i,2] - den[i,2])*x[2*i-1] + (num[i,1] - den[i,1])*x[2*i])/den[i,2] + uu[i]) else k[i]*((num[i,2] - den[i,2])*x[2*i-1] + (num[i,1] - den[i,1])*x[2*i] + uu[i]);
           elseif 2*(i - n_num2) <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
              der(x[2*i]) = if abs(den[i, 2])>eps and abs(num[i, 2])>eps then den[i, 2]*uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i] else uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i];
              uu[i+1] = if abs(den[i, 2])>eps and abs(num[2*i-n_num2-1, 2])>eps then k[i]*(((num[2*i-n_num2-1, 1]*num[2*i-n_num2-1 + 1, 1] - den[i, 2])*x[2*i-1] + (num[2*i-n_num2-1, 1] + num[2*i-n_num2-1 + 1, 1] - den[i, 1])*x[2*i])/den[i,2] + uu[i]) else k[i]*((num[2*i-n_num2-1, 1]*num[2*i-n_num2-1 + 1, 1] - den[i, 2])*x[2*i-1] + (num[2*i-n_num2-1, 1] + num[2*i-n_num2-1 + 1, 1] - den[i, 1])*x[2*i] + uu[i]);
           elseif 2*(i-n_num2) -1== n_num1 then
              // State space system in form (2) with 1 first order numerator polynomial
             der(x[2*i]) = if abs(den[i, 2])>eps then den[i, 2]*uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i] else uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i];
              uu[i+1] = if abs(den[i, 2])>eps then k[i]*(num[2*i-n_num2-1,1]*x[2*i-1] + x[2*i])/den[i,2] else k[i]*num[2*i-n_num2-1,1]*x[2*i-1] + x[2*i];
           else
              // State space system in form (3)
              der(x[2*i]) = if abs(den[i, 2])>eps then den[i, 2]*uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i] else uu[i] - den[i,2]*x[2*i-1] - den[i,1]*x[2*i];
              uu[i+1] = if abs(den[i, 2])>eps then k[i]*x[2*i-1]/den[i,2] else k[i]*x[2*i-1];
           end if;
        else
           // State space system in form (1) with 2 first order denominator polynomials
           der(x[2*i]) = if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then  den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]*uu[i] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])*x[2*i-1]  - (den[max(2*(i-n_den2)-1,i), 1]+den[max(2*(i-n_den2),i), 1])*x[2*i] else uu[i] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])*x[2*i-1]  - (den[max(2*(i-n_den2)-1,i), 1]+den[max(2*(i-n_den2),i), 1])*x[2*i];
           uu[i+1]     = if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>eps then k[i]*(((num[i,2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]))*x[2*i-1] + (num[i,1] - (den[max(2*(i-n_den2)-1,i), 1]+den[max(2*(i-n_den2),i), 1]))*x[2*i])/den[max(2*(i-n_den2)-1,i), 1]/den[max(2*(i-n_den2),i), 1] + uu[i]) else k[i]*(num[max(2*(i-n_num2),i), 2]*x[2*i-1] + (num[i,1] - (den[max(2*(i-n_den2),i), 1]))*x[2*i] + uu[i]);
        end if;
     end for;

     for i in i_d:n_den1 loop
        // Construct state space systems of order 1

        if n_num2 <= n_den2 and 2*(n_den2-n_num2)+i <= n_num1 then
           // State space system in form (4)
           der(x[2*n_den2+i]) =  if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1]*uu[i_k+i]-den[n_den2+i,1]*x[2*n_den2+i] else (num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1])*uu[i_k+i];
           uu[i_k+i+1] = if abs(den[n_den2 + i, 1])>eps then  k[i_k+i]*((num[max(1,n_num2 + 2*(n_den2-n_num2)+i),1]-den[n_den2+i,1])*x[2*n_den2+i]/den[n_den2+i,1] + uu[i_k+i]) else x[2*n_den2+i] + k[i_k+i]*uu[i_k+i];
        elseif n_num2 > n_den2 and i-i_d+1 <= n_num1 then
           // State space system in form (4)
           der(x[2*n_den2+i]) = if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1]*uu[i_k+i]-den[n_den2+i,1]*x[2*n_den2+i] else num[max(1, n_num2 + i - i_d + 1), 1]*uu[i_k+i];
           uu[i_k+i+1] = if abs(den[n_den2 + i, 1])>eps then k[i_k+i]*((num[max(1,n_num2 + i-i_d+1),1]-den[n_den2+i,1])*x[2*n_den2+i]/den[n_den2+i,1] + uu[i_k+i]) else x[2*n_den2+i] + k[i_k+i]*uu[i_k+i];
        else
           // State space system in form (5)
           der(x[2*n_den2+i]) = if abs(den[n_den2 + i, 1])>eps then den[n_den2 + i, 1]*uu[i_k+i]-den[n_den2+i,1]*x[2*n_den2+i] else uu[i_k+i];
           uu[i_k+i+1] = if abs(den[n_den2 + i, 1])>eps then k[i_k+i]*x[2*n_den2+i]/den[n_den2+i,1] else k[i_k+i]*x[2*n_den2+i];
        end if;
     end for;
     y = k_total*uu[i_k+n_den1+1];
 else
    for i in 1:size(uu, 1) - 1 loop
      uu[i + 1] = u;
    end for;
  end if;

 uu[1] = u;
  connect(x, discretePart.x);
  connect(y, discretePart.y[1]);

initial equation
  if continuous then
    if init ==Modelica_LinearSystems2.Controller.Types.Init.InitialState then
        x = x_start;
    elseif init ==Modelica_LinearSystems2.Controller.Types.Init.SteadyState then
        der(x) = zeros(nx);
    elseif init ==Modelica_LinearSystems2.Controller.Types.Init.InitialOutput and nx>0 then
        y = y_start;
        der(x[1:nx-(if nx>1 then 2 else 1)]) = zeros(nx-(if nx>1 then 2 else 1));

     end if;
  end if;
  annotation (
  defaultComponentName="zerosAndPoles",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Text(
          extent={{-88,-8},{90,52}},
          lineColor={0,0,127},
          textString="k*z(s)"),
        Line(points={{-82,-20},{78,-20}}, color={0,0,127}),
        Text(
          extent={{-92,-24},{88,-86}},
          lineColor={0,0,127},
          textString="p(s)"),
        Text(
          extent={{-100,96},{96,62}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-100,0},{-60,0}}, color={0,0,255}),
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,255}),
        Line(points={{60,0},{100,0}}, color={0,0,255}),
        Text(
          extent={{-52,50},{50,10}},
          lineColor={0,0,0},
          textString="k*z(s)"),
        Line(points={{50,0},{-50,0}}, color={0,0,0}),
        Text(
          extent={{-50,-10},{52,-52}},
          lineColor={0,0,0},
          textString="p(s)")}),
    Documentation(info="<html>
<p>
This function transforms a zeros-poles-gain system representation into state space representation.
To achieve well numerical condition the ZerosAndPoles transfer function is transformed into state space
form by creating first and second order blocks that are connected
together in series. Every block is represented in controller
canonical form and scaled such that the gain from the input
of this block to its output is one (i.e. y(s=0) = u(s=0)),
if this is possible. Details are given below.
</p>

<h4>Algorithmic details</h4>
<p>
The ZerosAndPoles transfer function is defined as:
</p>
<pre>         product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
  y = k*--------------------------------------------------------- * u
         product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre>
<p>
This is treated as a series connection of first and second order
systems. If size(n1) == size(d1) and size(n2) == size(d2)
this gives the following sequence of operations:
</p>
<pre>        s^2 + n2[1,1]*s + n2[1,2]
  y_1 = ------------------------- * u
        s^2 + d2[1,1]*s + d2[1,2]
&nbsp;
        s^2 + n2[2,1]*s + n2[2,2]
  y_2 = ------------------------- * y_1
        s^2 + d2[2,1]*s + d2[2,2]
&nbsp;
     ...
&nbsp;
        s + n1[..]
  y_n = ---------- * y_(n-1)
        s + d1[..]
&nbsp;
    y = k*y_n
</pre>
<p>
Based on this representation, evrey block with transfer function G(s) could be transformed into
</p>
<pre>  G(s) = k * F(s)
</pre>
<p>
with F(s) has unit gain. This leads to representations of the forms
</p>
<pre>           a2 + a1*s + s^2       a2      b2 + a1*b2/a2*s + b2/a2*s^2
  G(s) = -------------------- = ---- * ------------------------------ = k * F(s),  k = a2/b2  (1)
           b2 + b1*s + s^2       b2           b2 + b1*s + s^2
&nbsp;
for second order systems and
&nbsp;
           a + s     a     b + b/a*s
  G(s) = -------- = --- * ---------- = k * F(s),   k = a/b
           b + s     b      b + s
</pre>
<p>
for first order systems respectively.
</p>
<p>
The complete system is now considered as the series connections of all the single unit gain transfer functions and an overall gain k with
</p>
<pre>  k = product(ki).
</pre>
<p>
In the general case, the following system structures
and the corresponding state space systems can appear
(note, 'c' is the reciprocal local gain 1/k):
</p>
<pre>(1)
          a2 + a1*s + s^2           der(x1) = x2
    y = ---------------------  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
          b2 + b1*s + s^2                 y = c*((a2-b2)*x1 + (a1-b1)*x2 + u),  c = b2/a2
&nbsp;
(2)
             s + a                 der(x1) = x2
    y = ---------------- * u  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*s + s^2                  y = k*(a1/b2*x1 +x2/b2),  c = b2/a
&nbsp;
(3)
               1                  der(x1) = x2
    y = --------------- *u   -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*s + s^2                 y = c*x1/b2,  c = b2
&nbsp;
(4)
       a + s                       der(x) = -b*x + b*u
   y = ----- * u             -->        y = c*((a-b)/b*x + u),  c = b/a
       b + s
&nbsp;
(5)
         1
   y = ----- * u             -->   der(x) = -b*x + b*u
       b + s                            y = x,  c = b
</pre>
 <p>
If the sizes of the numerator and denominator polynomials
do not match, the small systems are built in the
following way:
</p>
<pre>(1) Build systems of form (1) by combining
    - 1 d2 and 1 n2
      (= 1 second order denominator and 1 second order numerator) or
    - 1 d2 and 2 n1 or
    - 2 d1 and 1 n2
(2) Build at most one system of form (2) by combining
    - 1 d2 and 1 n2
(3) Build systems of form (3) by
    - 1 d2
(4) Build systems of form (4) by combining
    - 1 d1 and 1 n1
(5) Build systems of form (5) by
    - 1 d1
</pre>
<p>
The numeric properties of the resulting state space system
depends on which first and second order polynomials are
combined and connected together. From a numerical point of view, it
would therefore be useful to combine the polynomials
based on the numeric values of the polynomial coefficients,
(e.g., in a first step the polynomials could be sorted
according to their cut-off frequency).
</p>
<p>
However, this has the disadvantage that the structure of the
resulting state space system depends on the numeric
values of the polynomial coefficients. Since Modelica
environments perform symbolic pre-processing on equations,
this would mean that a change of a polynomial coefficient
requires to newly compile the state space system.
</p>
<p>
If, on the other hand, the structure of the state
space system depends only on dimension information
of the n1,n2,d1,d2 arrays, then the polynomial coefficients
can be changed without a new translation of the model.
This is the major reason why the structure of the
state space system in the implementation of this block
is based only on dimension information.
</p>
<p>
This is, e.g., not critical for the provided filters:
The dimension of the n1,n2,d1,d2 arrays depend for
filters only on the filter characteristics
(Bessel, Butterworth etc.), the filter type (low pass,
high pass etc.) and on the filter order. If any
of this data is changed, the model has to be
newly compiled. All the other filter data, such as
cut-off frequency or ripple amplitude, can be changed
without re-compilation of the model.
The ZerosAndPoles transfer function is now constructed
for the filters in such a way that the filter zeros
and poles are appropriately sorted to give better
numerical properties.
</p>
<p>
Another alternative implementation of the state
space system would be to use the function controller canonical
form that directly results from the transfer function.
The severe disadvantage
of this approach is that the structure of the state
space system from above is lost for the symbolic preprocessing.
If, e.g., index reduction has to be applied (e.g. since a
filter is used to realize a non-linear inverse model),
then the tool cannot perform the index reduction.
Example:
</p>
<p>
Assume, a generic first order state space system
is present
</p>
<pre>   <b>der</b>(x) = a*x + b*u
        y = c*x + d*u
</pre>
<p>
and the values of the scalars a,b,c,d are parameters
that might be changed before the simulation starts.
If y has to be differentiated symbolically during code
generation, then
</p>
<pre>      <b>der</b>(y) = c*<b>der</b>(x) + d*<b>der</b>(u)
      <b>der</b>(x) = a*x + b*u
</pre>
<p>
As a result, u needs to be differentiated too, and this
might not be possible and therefore translation might fail.
</p>
<p>
On the other hand, if the first order system is
defined to be a low pass filter and the state space
system is generated by keeping this structure, we have
(see form (5) above):
</p>
<pre>   <b>der</b>(x) = -b*x + u
        y = x
</pre>
<p>
Differentiating y symbolically leads to:
</p>
<pre>     <b>der</b>(y) = <b>der</b>(x)
     <b>der</b>(x) = -b*x + u
</pre>
<p>
Therefore, in this case, the derivative of u is not
needed and the tool can continue with the symbolic
processing.
</p>
</html>"));
end ZerosAndPoles;
