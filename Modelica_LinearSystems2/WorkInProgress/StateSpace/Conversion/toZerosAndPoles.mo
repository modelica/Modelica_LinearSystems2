within Modelica_LinearSystems2.WorkInProgress.StateSpace.Conversion;
function toZerosAndPoles
  "Generate a zeros-and-poles representation from a SISO state space representation"

  import Modelica;
  import Modelica_LinearSystems2;
  import Complex;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.WorkInProgress;

  input StateSpace ss "StateSpace object";
  input Complex frequency2 = Complex(1);
//protected
//  input Boolean cancel=true "false to hinder cancellation";// is not fully realized
public
  output ZerosAndPoles zp;

protected
  StateSpace ssm= if size(ss.A,1)>0 then WorkInProgress.StateSpace.Transformation.toIrreducibleForm(ss) else StateSpace(ss.D[1,1]);
  Complex poles[:];
  Complex zeros[:];

  Real gain;

  Complex frequency;
  Complex Gs;
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

   zeros := WorkInProgress.StateSpace.Analysis.invariantZeros(ssm);

    if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
      assert(size(ss.B, 2) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of inputs is "
         + String(size(ss.B, 2)) + " instead of 1");
      assert(size(ss.C, 1) == 1, " function fromStateSpaceSISO expects a SISO-system as input\n but the number of outputs is "
         + String(size(ss.C, 1)) + " instead of 1");
    end if;
    zp := ZerosAndPoles(
        z=zeros,
        p=poles,
        k=1);

//    v := sum(cat(1, zeros[:].re,  -poles[:].re))/max(size(zeros,1)+size(poles,1),1)+13/23;
    v := sum(cat(1, zeros[:].re,  -poles[:].re))/max(size(poles,1),1)+13/23;
    frequency := Complex(v)*19/17;
    Modelica.Utilities.Streams.print(String(v));

    Gs := ZerosAndPoles.Analysis.evaluate(zp, frequency);

    As := -ssm.A;
    for i in 1:size(As, 1) loop
      As[i, i] := As[i, i] + frequency.re;
    end for;

    pk := StateSpace.Internal.partialGain(As, ssm.B[:, 1]);
    gain := (ssm.C[1, size(As, 1)]*pk + ss.D[1, 1])/Gs.re;

    zp := ZerosAndPoles(
        z=zeros,
        p=poles,
        k=gain);

  else
    zp := ZerosAndPoles(
        z=fill(Complex(0), 0),
        p=fill(Complex(0), 0),
        k=scalar(ss.D));

  end if;
  zp.uName := ss.uNames[1];
  zp.yName := ss.yNames[1];

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
zp = StateSpace.Conversion.<b>toZerosAndPoles</b>(ss)
</pre> </blockquote>

<h4>Description</h4>
<p>
Computes a ZerosAndPoles record
</p>
<blockquote><pre>
          product(s + n1[i]) * product(s^2 + n2[i,1]*s + n2[i,2])
zp = k * ---------------------------------------------------------
          product(s + d1[i]) * product(s^2 + d2[i,1]*s + d2[i,2])
</pre></blockquote>
<p>
of a system from state space representation using the transformation algorithm described in [1].
The uncontrollable and unobservable parts are isolated and the eigenvalues and invariant zeros of the controllable and observable sub system are calculated.
</p>

<h4>Example</h4>
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

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga, A and Sima, V. (1981):</dt>
<dd> <b>Numerically stable algorithm for transfer function matrix evaluation</b>.
     Int. J. Control, Vol. 33, No. 6, pp. 1123-1133.<br>&nbsp;</dd>
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
end toZerosAndPoles;
