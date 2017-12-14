within Modelica_LinearSystems2.WorkInProgress.ZerosAndPoles.Conversion;
function toStateSpace2
  "Transform a ZerosAndPoles object into a StateSpace object"
 //encapsulated function fromZerosAndPoles
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.Math.Vectors;
  import Modelica_LinearSystems2.StateSpace;
  import Complex;

  input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
  output StateSpace ss(
    redeclare Real A[ZerosAndPoles.Analysis.denominatorDegree(zp),
      ZerosAndPoles.Analysis.denominatorDegree(zp)],
    redeclare Real B[ZerosAndPoles.Analysis.denominatorDegree(zp),1],
    redeclare Real C[1,ZerosAndPoles.Analysis.denominatorDegree(zp)],
    redeclare Real D[1,1]) "Transfer function in StateSpace SISO form";

protected
  Real A[2,2] "system matrix of partial 2nd order system";
  Real B[2,1] "input matrix of partial 2nd order system";
  Real C[1,2] "output matrix of partial 2nd order system";
  Real D[1,1] "feedthrough matrix of partial 2nd order system";
  Real a "system 'matrix' of partial 1st order system";
  Real b "input 'matrix' of partial 1st order system";
  Real c "output 'matrix' of partial 1st order system";
  Real d "feedthrough 'matrix' of partial 1st order system";
  Integer nx=max(ZerosAndPoles.Analysis.numeratorDegree(zp),ZerosAndPoles.Analysis.denominatorDegree(zp));
  Integer n_num1=size(zp.n1, 1);
  Integer n_num2=size(zp.n2, 1);
  Integer n_den1=size(zp.d1, 1);
  Integer n_den2=size(zp.d2, 1);
  Integer n_num=n_num1 + 2*n_num2;
  Integer n_den=n_den1 + 2*n_den2;

  Integer i_d=if n_num2 > n_den2 then 2*(n_num2 - n_den2) + 1 else 1;
  Integer i_k=if n_num2 > n_den2 then n_den2 - (n_num2 - n_den2) else n_den2;
  Integer i;
  Integer ili=if max(n_den2, n_num2) > 0 then i_d else  max(2, i_d);

  Real num[nx,2]=[zp.n2; [zp.n1,zeros(n_num1)]; zeros(max(0,nx - n_num2 - n_num1), 2)]
    "Numerator matrix, in order that indices are defined in all situations in all if clauses";
  Real den[nx,2]=[zp.d2; [zp.d1,zeros(n_den1)]; zeros(max(0,nx - n_den2 - n_den1), 2)]
    "Denominator matrix, in order that indices are defined in all situations in all if clauses";
  Real k[i_k + n_den1](each fixed=false)
    "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
  Real k_total;

  Boolean dZero=true;

algorithm
  assert(n_num <= n_den,
    "ZerosAndPoles transfer function is not proper as required from StateSpace system:\n"
     + "  numerator degree (= " + String(n_num) +
    ") <= denominator degree (= " + String(n_den) + ") required.");

  if n_den > 0 then
    for i in 1:max(n_den2, n_num2) loop
      // State space systems of order 2
      if i <= n_den2 then
        if i <= n_num2 then
            // State space system in form (1)
          k[i] := StateSpace.Internal.scaleFactor2(
              num[i, 1],
              num[i, 2],
              den[i, 1],
              den[i, 2]);
        elseif 2*(i - n_num2) <= n_num1 then
            // State space system in form (1) with 2 first order numerator polynomials
          k[i] := StateSpace.Internal.scaleFactor2(
              num[2*(i - n_num2)-1, 1] + num[2*(i - n_num2), 1],
              num[2*(i - n_num2)-1, 1]*num[2*(i - n_num2), 1],
              den[i, 1],
              den[i, 2]);
        elseif  2*(i-n_num2) -1== n_num1 then
            // State space system in form (2) with 1 first order numerator polynomial
//          k[i] := StateSpace.Internal.scaleFactor2(
//              num[2*i-n_num2-1, 1],
//              0,
//              den[i, 1],
//              den[i, 2]);
          k[i] := StateSpace.Internal.scaleFactor2(
              1,
              num[2*i-n_num2-1, 1],
              den[i, 1],
              den[i, 2]);
        else
            // State space system in form (3)
          k[i] := StateSpace.Internal.scaleFactor2(
              1,
              1,
              den[i, 1],
              den[i, 2]);
        end if;
      else
         // State space system in form (1) with 2 first order denominator polynomials
        k[i] := StateSpace.Internal.scaleFactor2(
            num[i, 1],
            num[i, 2],
            den[2*(i - n_den2)-1, 1] + den[2*(i - n_den2), 1],
            den[2*(i - n_den2)-1, 1]*den[2*(i - n_den2), 1]);
      end if;
    end for;

    for i in i_d:n_den1 loop
      // State space systems of order 1
      if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
         // State space system in form (4)
        k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
          n_num2) + i), 1], den[n_den2 + i, 1]);
      elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
         // State space system in form (4)
        k[i_k + i] := StateSpace.Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
          1], den[n_den2 + i, 1]);
      else
         // State space system in form (5)
        k[i_k + i] := StateSpace.Internal.scaleFactor1(1, den[n_den2 + i, 1]);
      end if;
    end for;

    k_total := zp.k/product(k);

    ss.A := zeros(nx, nx);
    ss.B := zeros(nx, 1);
    ss.C := zeros(1, nx);
    ss.D := zeros(1, 1);

 // Calculation of matrices A, B, C, D
 //first elements of A, B, C and D
    if max(n_den2, n_num2) > 0 then
      A[1, :] := {0,1};
      B[1, 1] := 0;
      // Construct state space systems of order 2
      if 1 <= n_den2 then
        A[2, :] := {-den[1, 2],-den[1, 1]};
        B[2, 1] := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then den[1, 2] else 1;

        if 1 <= n_num2 then
         // State space system in form (1)
          C := if abs(den[1, 2])>Modelica.Constants.eps and abs(num[1, 2])>Modelica.Constants.eps then [num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/num[1, 2] else [num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
          D := [1/k[1]];
          dZero := false;
       elseif 1 - n_num2 + 1 <= n_num1 then
        // State space system in form (1) with 2 first order numerator polynomials
        B[2, 1] := if abs(den[1, 2])>Modelica.Constants.eps then den[1, 2] else 1;
          C := if abs(den[1, 2])>Modelica.Constants.eps then   [num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/num[1, 1]/num[1, 2] else [num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
          D := [1/k[1]];
          dZero := false;
       elseif 1 - n_num2 == n_num1 then
        // State space system in form (2) with 1 first order numerator polynomial
          B[2, 1] := if abs(den[1, 2])>Modelica.Constants.eps then den[1, 2] else 1;
          C := if abs(num[1, 1])>Modelica.Constants.eps then [1,1/num[1, 1]] else [num[1, 1],1];
          D := [0];
          dZero := dZero and true;
        else
               // State space system in form (3)
          B[2, 1] := if abs(den[1, 2])>Modelica.Constants.eps then den[1, 2] else 1;
          C := if abs(den[1, 2])>Modelica.Constants.eps then  [1,0]/num[1, 2] else [1,0];
          D := [0];
          dZero := dZero and true;
        end if;
      else
        // State space system in form (1) with 2 first order denominator polynomials
        A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
        B[2, 1] := if abs(den[1, 1]*den[2, 1])>Modelica.Constants.eps then den[1, 1]*den[2, 1] else 1;
        C := if abs(den[1, 1]*den[2, 1])>Modelica.Constants.eps then  [num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/num[1, 1]/num[2, 1] else [num[1, 2],num[1, 1] - den[2, 1]];
        D := [1/k[1]];
        dZero := false;
      end if;
      ss.A[1:2, 1:2] := A;
      ss.B[1:2, 1] := vector(B);
      ss.C[1, 1:2] := vector(C);
      ss.D := D;

    else

   // Construct state space systems of order 1
      a := -den[1, 1];

      if 1 <= n_num1 then
        // State space system in form (4)
        b := if abs(den[1, 1])>Modelica.Constants.eps then den[1,1] else num[1,1];
        c := if abs(den[1, 1])>Modelica.Constants.eps then k[1]*(num[1, 1] - den[1, 1])/den[1, 1] else k[1];
        d := k[1];
        dZero := false;
      else
       // State space system in form (5)
        b := if abs(den[1, 1])>Modelica.Constants.eps then den[1,1] else if n_num1>0 then num[1,1] else 1;
        c := if abs(den[1, 1])>Modelica.Constants.small then 1/num[1, 1] else 1;
        d := 0;
        dZero := dZero and true;
      end if;
      ss.A[1, 1] := a;
      ss.B[1, 1] := b;
      ss.C[1, 1] := c;
      ss.D[1, 1] := d;

    end if;
 /// for i=2 to degree(system)
    A[1, :] := {0,1};
    B[1, 1] := 0;
    for i in 2:max(n_den2, n_num2) loop
         // Construct state space systems of order 2
      if i <= n_den2 then
        A[2, :] := {-den[i, 2],-den[i, 1]};
        B[2, 1] := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then den[i, 2] else 1;

        if i <= n_num2 then
               // State space system in form (1)

          C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[i, 2])>Modelica.Constants.eps then [num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/num[i, 2] else [num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
          D := [1/k[i]];
          dZero := false;

        elseif 2*(i - n_num2) <= n_num1 then
          // State space system in form (1) with 2 first order numerator polynomials
          C := if abs(den[i, 2])>Modelica.Constants.eps and abs(num[2*i-n_num2-1, 2])>Modelica.Constants.eps then [num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]]/num[i, 2] else [num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]];
          D := [1/k[i]];
          dZero := false;

        elseif 2*(i-n_num2) -1== n_num1 then
        // State space system in form (2) with 1 first order numerator polynomial
          B[2, 1] := if abs(den[i, 2])>Modelica.Constants.eps then den[i, 2] else 1;
          C := if abs(den[i, 2])>Modelica.Constants.eps then [num[2*i-n_num2-1, 1],1]/num[i, 2] else [num[2*i-n_num2-1, 1],1];
          D := [0];
          dZero := dZero and true;
        else
          // State space system in form (3)
          B[2, 1] := if abs(den[i, 2])>Modelica.Constants.eps then den[i, 2] else 1;
          C := if abs(den[i, 2])>Modelica.Constants.eps then [1,0]/num[i, 2] else [1,0];
          D := [0];
          dZero := dZero and true;
        end if;

      else
            // State space system in form (1) with 2 first order denominator polynomials
        A[2, :] := {-(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),-(den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])};
        B[2, 1] := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>Modelica.Constants.eps then den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1] else 1;
        C := if abs(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1])>Modelica.Constants.eps then [num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/num[max(2*(i-n_den2)-1,i), 1]/num[max(2*(i-n_den2),i), 1] else [num[max(2*(i-n_num2),i), 2],num[max(2*(i-n_num2),i), 1] - den[max(2*(i-n_den2),i), 1]];
        D := [1/k[i]];
        dZero := false;
      end if;
      ss.A[2*i, 1:2*i - 2] := B[2, 1]*ss.C[1, 1:2*i - 2];
      ss.A[2*i - 1:2*i, 2*i - 1:2*i] := A;
      ss.B[2*i, 1] := if dZero then 0 else B[2, 1]*ss.D[1, 1];
      ss.C[1, 1:2*i - 2] := if dZero then fill(0, 2*i - 2) else D[1, 1]*ss.C[
        1, 1:2*i - 2];
      ss.C[1, 2*i - 1:2*i] := vector(C);
      ss.D := D*ss.D;
    end for;
 //  for i in max(2,i_d):n_den1 loop
    for i in ili:n_den1 loop
         // Construct state space systems of order 1
      a := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then -den[n_den2 + i, 1] else 0.0;

      if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
            // State space system in form (4)

        c := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then (num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/num[n_den2 + i, 1] else 1.0;
        b := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then den[n_den2 + i, 1] else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1];
        d := 1/k[i_k + i];
        dZero := false;
      elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
      // State space system in form (4)
        c := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then (num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/num[n_den2 + i, 1] else 1.0;
        b := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then den[n_den2 + i, 1] else num[max(1, n_num2 + i - i_d + 1), 1];
        d := 1/k[i_k + i];
        dZero := false;
      else
       // State space system in form (5)
        c := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then 1/num[n_den2 + i, 1] else 1;
        b := if abs(den[n_den2 + i, 1])>Modelica.Constants.eps then den[n_den2 + i, 1] else 1;
        d := 0;
        dZero := dZero and true;
      end if;
      ss.A[2*n_den2 + i, 1:2*n_den2 + i - 1] := b*ss.C[1, 1:2*n_den2 + i - 1];
      ss.A[2*n_den2 + i, 2*n_den2 + i] := a;
      ss.B[2*n_den2 + i, 1] := if dZero then 0 else b*ss.D[1, 1];
      ss.C[1, 1:2*n_den2 + i - 1] := if dZero then fill(0, 2*n_den2 + i - 1) else
              d*ss.C[1, 1:2*n_den2 + i - 1];
      ss.C[1, 2*n_den2 + i] := c;
      ss.D := if dZero then [0] else d*ss.D;

    end for;

    ss.C := k_total*ss.C;
    ss.D := k_total*ss.D;
  else
    ss := Modelica_LinearSystems2.StateSpace(zp.k);
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
ss = ZerosAndPoles.Conversion.toStateSpace<b>toStateSpace</b>(zp)
</pre></blockquote>

<h4>Description</h4>
<p>
This function transforms a zeros-poles-gain system representation into state space representation.
To achieve well numerical condition the ZerosAndPoles transfer function is transformed into state space
form by creating first and second order blocks that are connected
together in series. Every block is represented in controller
canonical form and scaled such that the gain from the input
of this block to its output is one (i.e. y(p=0) = u(p=0)),
if this is possible. Details are given below.
</p>

<h4>Algorithmic details</h4>
<p>
The ZerosAndPoles transfer function is defined as:
</p>
<blockquote><pre>
         product(p + n1[i]) * product(p^2 + n2[i,1]*p + n2[i,2])
y = k * --------------------------------------------------------- * u
         product(p + d1[i]) * product(p^2 + d2[i,1]*p + d2[i,2])
</pre></blockquote>
<p>
This is treated as a series connection of first and second order
systems. If size(n1) == size(d1) and size(n2) == size(d2)
this gives the following sequence of operations:
</p>
<blockquote><pre>
      p^2 + n2[1,1]*p + n2[1,2]
y_1 = ------------------------- * u
      p^2 + d2[1,1]*p + d2[1,2]
&nbsp;
      p^2 + n2[2,1]*p + n2[2,2]
y_2 = ------------------------- * y_1
      p^2 + d2[2,1]*p + d2[2,2]
&nbsp;
  ...
&nbsp;
      p + n1[..]
y_n = ---------- * y_(n-1)
      p + d1[..]
&nbsp;
  y = k*y_n
</pre></blockquote>
<p>
Based on this representation, evrey block with transfer function G(p) could be transformed into
</p>
<blockquote><pre>
G(p) = k * F(p)
</pre></blockquote>
<p>
with F(p) has unit gain. This leads to representations of the forms
</p>
<blockquote><pre>
        a2 + a1*p + p^2     a2     b2 + a1*b2/a2*p + b2/a2*p^2
G(p) = ----------------- = ---- * ----------------------------- = k * F(p),  k = a2/b2  (1)
        b2 + b1*p + p^2     b2           b2 + b1*p + p^2
</pre></blockquote>
<p>
for second order systems and
</p>
<blockquote><pre>
        a + p     a     b + b/a*p
G(p) = ------- = --- * ----------- = k * F(p),   k = a/b
        b + p     b      b + p
</pre></blockquote>
<p>
for first order systems respectively.
</p>
<p>
The complete system is now considered as the series connections of all
the single unit gain transfer functions and an overall gain k with
</p>
<blockquote><pre>
k = product(ki).
</pre></blockquote>
<p>
In the general case, the following system structures
and the corresponding state space systems can appear
(note, 'c' is the reciprocal local gain 1/k):
</p>
<blockquote><pre>
(1)
          a2 + a1*p + p^2           der(x1) = x2
    y = ---------------------  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
          b2 + b1*p + p^2                 y = c*((a2-b2)/b2*x1 + (a1-b1)/b2*x2 + u),  c = b2/a2
&nbsp;
(2)
             p + a                 der(x1) = x2
    y = ---------------- * u  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                  y = k*(a/b2*x1 +x2/b2),  c = b2/a
&nbsp;
(3)
               1                  der(x1) = x2
    y = --------------- *u   -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                 y = c*x1/b2,  c = b2
&nbsp;
(4)
       a + p                       der(x) = -b*x + b*u
   y = ----- * u             -->        y = c*((a-b)/b*x + u),  c = b/a
       b + p
&nbsp;
(5)
         1
   y = ----- * u             -->   der(x) = -b*x + b*u
       b + p                            y = x,  c = b

</pre></blockquote>

<p>
If the sizes of the numerator and denominator polynomials
do not match, the small systems are built in the
following way:
</p>
<blockquote><pre>
(1) Build systems of form (1) by combining
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
</pre></blockquote>
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
<blockquote><pre>
<b>der</b>(x) = a*x + b*u
     y = c*x + d*u
</pre></blockquote>
<p>
and the values of the scalars a,b,c,d are parameters
that might be changed before the simulation starts.
If y has to be differentiated symbolically during code
generation, then
</p>
<blockquote><pre>
<b>der</b>(y) = c*<b>der</b>(x) + d*<b>der</b>(u)
<b>der</b>(x) = a*x + b*u
</pre></blockquote>
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
<blockquote><pre>
<b>der</b>(x) = -b*x + u
      y = x
</pre></blockquote>
<p>
Differentiating y symbolically leads to:
</p>
<blockquote><pre>
<b>der</b>(y) = <b>der</b>(x)
<b>der</b>(x) = -b*x + u
</pre></blockquote>
<p>
Therefore, in this case, the derivative of u is not
needed and the tool can continue with the symbolic
processing.
</p>

<h4>Example</h4>
<blockquote><pre>
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2 + p +1);

<b>algorithm</b>
  ss := Modelica_LinearSystems2.ZerosAndPoles.Conversion.toStateSpace(zp);
// ss.A = [0, 1; -1, -1],
// ss.B = [0; 1],
// ss.C = [1, 1],
// ss.D = [0],
</pre></blockquote>
</html>"));
end toStateSpace2;
