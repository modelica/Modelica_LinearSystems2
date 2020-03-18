within Modelica_LinearSystems2.Internal;
encapsulated function frequencyEvaluate
  "Evaluate a SISO transfer function defined by Zeros and Poles matrices at a given complex value re + j*im"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Internal;
  import Modelica.Utilities.Streams.print;

  input Real gain "Gain of transfer function";
  input Real Zeros[:,2]
    "Zeros as Real matrix (first column: real, second column imaginary values)";
  input Real Poles[:,2]
    "Poles as Real matrix (first column: real, second column imaginary values)";
  input Real re "Real value of complex number";
  input Real im "Imaginary value of complex number";
  output Real A "Amplitude of complex value";
  output Modelica.Units.SI.Angle phi "Angle of complex value";
  output Integer info
    "= 0/1/2 success/infinity(A is a large value)/indefinite (A=0/0; A=0 returned)";
protected
  Integer nn = size(Zeros,1);
  Integer nd = size(Poles,1);
  Real Zeros2[size(Zeros,1)];
  Real n_A[ size(Zeros,1)];
  Real n_A2[size(Zeros,1)];
  Real d_A[ size(Poles,1)];
  Real d_A2[size(Poles,1)];
  Modelica.Units.SI.Angle n_phi[size(Zeros, 1)];
  Modelica.Units.SI.Angle d_phi[size(Poles, 1)];
  Boolean n_zero=false;
  Boolean d_zero=false;
algorithm
  assert( nd >= 1, "Vector Poles must have at least one element");
  assert( nn <= nd, "Size of vector Zeros <= size of vector Poles required");

  // Compute values of zeros and poles at s=re+j*im
  // z(s) = s - (a+j*b); z(re+j*im) = re-a + j(im-b)
  if nn > 0 then
     (n_A, n_phi, n_zero) :=Internal.toPolarForm([fill(re,nn)-Zeros[:, 1],fill(im,nn)-Zeros[:, 2]]);
  end if;
  (d_A, d_phi, d_zero) :=Internal.toPolarForm([fill(re,nd)-Poles[:, 1],fill(im,nd)-Poles[:, 2]]);

  // Handle zeros in numerator
  if n_zero then
     A := 0;
     phi := 0;
     info :=if d_zero then 2 else 0;
     return;
  end if;

   // Compute angle of fraction
   phi :=0;
   for i in 1:nn loop
      phi := phi + n_phi[i] - d_phi[i];
   end for;

   for i in nn+1:nd loop
      phi := phi - d_phi[i];
   end for;

   // Compute absolute value (avoid overflow)
   if d_zero then
      info :=1;
      A :=Modelica.Constants.inf;
   else
      info :=0;
      n_A2 :=Modelica.Math.Vectors.sort(n_A, ascending=false);
      d_A2 :=Modelica.Math.Vectors.sort(d_A, ascending=false);
      A := 1;
      for i in 1:nn loop
         A :=A*(n_A2[i]/d_A2[i]);
      end for;
      for i in nn+1:nd loop
         A :=A/d_A2[i];
      end for;
   end if;

  // Take gain into account
  if info == 0 then
     A :=gain*A;
  elseif info == 1 then
     A :=if gain >= 0 then A else -A;
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = ZerosAndPoles.Analysis.<b>evaluate</b>(zp, p, den_min=0)
</pre></blockquote>

<h4>Description</h4>
<p>
Function Analysis.<b>evaluate</b> evaluates the ZerosAndPoles transfer function at a given (complex) value of p and returns the value G(p)=N(p)/D(p). The optional argument den_min with default 0 is used to guard against a division by zero.
</p>
<pre>
  <b>if</b> |(D(p))| >= den_min <b>then</b>
     G(p) = N(p) / D(p);
  <b>elseif</b> D(p).re >= 0.0 <b>then</b>
     G(p) = N(p) / den_min
  <b>else</b>
     G(p) = -N(p) / den_min
  <b>end if</b>;
</p>
</pre>

<h4>Example</h4>
<blockquote><pre>
  Complex j = Modelica_LinearSystems2.Math.Complex.j();
  ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
  Modelica_LinearSystems2.ZerosAndPoles zp=(p+1)/(p^2+p+1);

  Complex result;

<b>algorithm</b>
  result := Modelica_LinearSystems2.ZerosAndPoles.Analysis.evaluate(zp, j+1);
//  result = 0.538462 - 0.307692j
</pre></blockquote>

<h4>See also</h4>
<p>
<a href=\"modelica://Modelica_LinearSystems2.Math.Polynomial.evaluateComplex\">Math.Polynomial.evaluateComplex</a>
</p>
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
end frequencyEvaluate;
