within Modelica_LinearSystems2.Internal;
function balanceABC "Return a balanced form of a system [A,B;C,0] to improve its condition by a state transformation (this function is provided in the trunk version
  of MSL and will be removed once Modelica_LinearSystems2 is based on the next MSL version)"
  extends Modelica.Icons.Function;

  input Real A[:, size(A, 1)] "System matrix A";
  input Real B[size(A,1),:] = fill(0.0, size(A,1),0)
    "System matrix B (need not be present)";
  input Real C[:,size(A,1)] = fill(0.0, 0, size(A,1))
    "System matrix C (need not be present)";
  output Real scale[size(A, 1)]
    "diagonal(scale)=T is such that [inv(T)*A*T, inv(T)*B; C*T, 0] has smaller condition as [A,B;C,0]";
  output Real As[size(A, 1), size(A, 1)] "Balanced matrix A (= inv(T)*A*T )";
  output Real Bs[size(A, 1), size(B, 2)] "Balanced matrix B (= inv(T)*B )";
  output Real Cs[size(C, 1), size(A, 1)] "Balanced matrix C (= C*T )";

protected
  Integer na=size(A, 1);
  Integer radix=2 "Radix of exponent representation must be 'radix'
          or a multiple of 'radix'";
  Integer radix2=radix*radix;
  Boolean noconv=true;
  Integer i=1;
  Integer j=1;
  Real CO;
  Real RO;
  Real G;
  Real F;
  Real S;
algorithm
  scale := ones(na);
  As := A;
  Bs := B;
  Cs := C;
  while noconv loop
    noconv := false;
    for i in 1:na loop
      CO := sum(abs(As[:, i])) - abs(As[i, i]) + sum(abs(Cs[:,i]));
      RO := sum(abs(As[i, :])) - abs(As[i, i]) + sum(abs(Bs[i,:]));
      G := RO/radix;
      F := 1;
      S := CO + RO;
      while not (CO >= G or CO == 0) loop
        F := F*radix;
        CO := CO*radix2;
      end while;
      G := RO*radix;
      while not (CO < G or RO == 0) loop
        F := F/radix;
        CO := CO/radix2;
      end while;
      if not ((CO + RO)/F >= 0.95*S) then
        G := 1/F;
        scale[i] := scale[i]*F;
        As[i, :] := As[i, :]*G;
        As[:, i] := As[:, i]*F;
        Bs[i, :] := Bs[i, :]*G;
        Cs[:, i] := Cs[:, i]*F;
        noconv := true;
      end if;
    end for;
  end while;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(scale,As,Bs,Cs) = Matrices.<strong>balanceABC</strong>(A,B,C);
(scale,As,Bs)    = Matrices.<strong>balanceABC</strong>(A,B);
(scale,As,,Cs)   = Matrices.<strong>balanceABC</strong>(A,C=C);
</pre></blockquote>

<h4>Description</h4>

<p>
This function returns a&nbsp;vector <strong>scale</strong>, such that with
<strong>T</strong>&nbsp;=&nbsp;diagonal(<strong>scale</strong>) system matrix <strong>S_scale</strong>
</p>
<blockquote>
  <table style=\"border:0\">
    <tr>
      <td>   </td>
      <td> | </td>
      <td style=\"text-align:center;\"> inv(<strong>T</strong>)*<strong>A</strong>*<strong>T</strong> </td>
      <td>   </td>
      <td style=\"text-align:center;\"> inv(<strong>T</strong>)*<strong>B</strong> </td>
      <td> | </td>
    </tr>
    <tr>
      <td> <strong>S_scale</strong> = </td>
      <td> | </td>
      <td>   </td>
      <td> &ensp;  </td>
      <td>   </td>
      <td> | </td>
    </tr>
    <tr>
      <td>   </td>
      <td> | </td>
      <td style=\"text-align:center;\"> <strong>C</strong>*<strong>T</strong> </td>
      <td>   </td>
      <td style=\"text-align:center;\"> <strong>0</strong> </td>
      <td> | </td>
    </tr>
  </table>
</blockquote>
<p>
has a better condition then system matrix&nbsp;<strong>S</strong>
</p>
<blockquote>
  <table style=\"border:0\">
    <tr>
      <td>   </td>
      <td> | <strong>A</strong> &ensp; </td>
      <td style=\"text-align:right;\"> <strong>B</strong> | </td>
    </tr>
    <tr>
      <td> <strong>S</strong> = </td>
      <td> | </td>
      <td style=\"text-align:right;\"> | </td>
    </tr>
    <tr>
      <td>   </td>
      <td> | <strong>C</strong> </td>
      <td style=\"text-align:right;\"> <strong>0</strong> |</td>
    </tr>
  </table>
</blockquote>
<p>
that is, conditionNumber(<strong>S_scale</strong>)&nbsp;&le;&nbsp;conditionNumber(<strong>S</strong>).
The elements of vector <strong>scale</strong> are multiples of&nbsp;2 which means that this function
does not introduce round-off errors.
</p>

<p>
Balancing a linear dynamic system in state space form
</p>

<blockquote><pre>
der(<strong>x</strong>) = <strong>A</strong>*<strong>x</strong> + <strong>B</strong>*<strong>u</strong>
    <strong>y</strong>  = <strong>C</strong>*<strong>x</strong> + <strong>D</strong>*<strong>u</strong>
</pre></blockquote>

<p>
means to find a state transformation
<strong>x_new</strong>&nbsp;= <strong>T</strong>*<strong>x</strong>&nbsp;= diagonal(<strong>scale</strong>)*<strong>x</strong>
so that the transformed system is better suited for numerical algorithms.
</p>

<h4>Example</h4>

<blockquote><pre>
<strong>import</strong> Modelica.Math.Matrices;

A = [1, -10,  1000; 0.01,  0,  10; 0.005,  -0.01,  10];
B = [100, 10; 1,0; -0.003, 1];
C = [-0.5, 1, 100];

(scale, As, Bs, Cs) := Matrices.balanceABC(A,B,C);
T    = diagonal(scale);
Diff = [Matrices.inv(T)*A*T, Matrices.inv(T)*B;
        C*T, zeros(1,2)] - [As, Bs; Cs, zeros(1,2)];
err  = Matrices.norm(Diff);

// scale = {16, 1, 0.0625}
// norm(A)  = 1000.15, norm(B)  = 100.504, norm(C)  = 100.006
// norm(As) = 10.8738, norm(Bs) = 16.0136, norm(Cs) = 10.2011
// err = 0
</pre></blockquote>

<p>
The algorithm is taken from
</p>
<dl>
<dt>H. D. Joos, G. Grbel:
<dd><strong>RAsP'91 Regulator Analysis and Synthesis Programs</strong><br>
    DLR - Control Systems Group 1991
</dl>
<p>
which is based on the <code>balance</code> function from EISPACK.
</p>
</html>",
        revisions="<html>
<ul>
<li><em>Sept. 14, 2014</em>
       by Martin Otter: Implemented.
</li>
</ul>
</html>"));
end balanceABC;
