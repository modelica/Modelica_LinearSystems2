within Modelica_LinearSystems2.WorkInProgress.ZerosAndPoles.Conversion;
function toTransferFunction
  "Generate a TransferFunction object from a ZerosAndPoles object"
  //encapsulated function fromZerosAndPoles
  import Modelica;
  import Complex;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input ZerosAndPoles zp "ZerosAndPoles transfer function of a system";
  output TransferFunction tf(redeclare Real n[2*size(zp.n2,1)+size(zp.n1,1)+1], redeclare Real
           d[                                                                                    2*size(zp.d2,1)+size(zp.d1,1)+1]);

protected
  Real k;
  Complex z[:];
  Complex p[:];
  Polynomial pn;
  Polynomial pd;
algorithm
  (z,p,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
  pn := Polynomial(z)*Polynomial(k);
  pd := Polynomial(p);
  tf.n := pn.c;
  tf.d := pd.c;
  tf.uName := zp.uName;
  tf.yName := zp.yName;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<table>
<tr> <td align=right>  tf </td><td align=center> =  </td>  <td> ZerosAndPoles.Conversion.<b>toTransferFunction</b>(zp)  </td> </tr>
</table>
<h4>Description</h4>
<p>
Computes a TransferFunction record
 <blockquote><pre>
           n(s)     b0 + b1*s + ... + bn*s^n
   tf = -------- = --------------------------
           d(s)     a0 + a1*s + ... + an*s^n
 </pre></blockquote>
from a ZerosAndPoles record representated by first and second order numerator and denominator polynomials. The poles and zeros and the gain <tt>k</tt> are computed (<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Analysis.zerosAndPoles\">zerosAndPoles</a>) and are used as inputs in the TransferFunction constructor.


<h4>Example</h4>
<blockquote><pre>
   ZerosAndPoles p = Modelica_LinearSystems2.ZerosAndPoles.p();
   Modelica_LinearSystems2.ZerosAndPoles zp = 1/(p + 3)/(p + 1)


<b>algorithm</b>
  tf:=Modelica_LinearSystems2.ZerosAndPoles.Conversion.toTransferFunction(zp);
//  tf = 1/( s^2 + 4*s + 3 )
</pre></blockquote>



</html>"));
end toTransferFunction;
