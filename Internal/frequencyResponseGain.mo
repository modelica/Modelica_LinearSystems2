within Modelica_LinearSystems2.Internal;
encapsulated function frequencyResponseGain
  "Compute the gain of a frequency response (system must be a SISO system)"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Internal;
  import Modelica.Utilities.Streams.print;

  input Real A[:,size(A,1)] "A-matrix of linear state space system";
  input Real B[size(A,1),:] "B-matrix of linear state space system";
  input Real C[:,size(A,1)] "C-matrix of linear state space system";
  input Real D[size(C,1),size(B,2)] "D-matrix of linear state space system";
  input Real Zeros[:,2]
    "Zeros of state space system as Real matrix (first column: real, second column imaginary values)";
  input Real Poles[:,2]
    "Poles of state space system as Real matrix (first column: real, second column imaginary values)";
  output Real gain
    "y(s) = gain*(s-z1)*(s-z2)*...*(s-zm)/((s-p1)*(s-p2)*...(s-pn))";

protected
  Real r;
  Real phi;
  Real kr1;
  Real kr2;
  Real A2[size(A,1), size(A,2)];

  function getReOutsidePolesZeros
    "Get p on the real axis so that there is a minimum distance to all poles and zeros"
    input Real Zeros[:,2]
      "Zeros of ss as Real matrix (first column: real, second column imaginary values)";
    input Real Poles[:,2]
      "Poles of ss as Real matrix (first column: real, second column imaginary values)";
    input Real minimumDistance=0.1;
    output Real p
      "Value on real axis > 0.0, so that poles and zeros on the axis have a minimumDistance to it";
    /* Most systems have no or only a few unstable poles or zeros.
    Searching for a suitable p is therefore fastest when searching 
    only in the unstable region, that is p > 0.0
    */
  protected
    Integer nVec=size(Poles, 1) + size(Zeros, 1);
    Real vec[:];
    Real vecSorted[:];
    Integer i;
    Integer iMax;
    Real small = minimumDistance*1e-6;
  algorithm
    i := 0;
    vec := zeros(nVec);
    for j in 1:size(Poles, 1) loop
      if Poles[j,1] > small and abs(Poles[j,2]) <= minimumDistance then
        i := i + 1;
        vec[i] := Poles[j,1];
      end if;
    end for;
    for j in 1:size(Zeros, 1) loop
      if Zeros[j,1] > small and abs(Zeros[j,2]) <= minimumDistance then
        i := i + 1;
        vec[i] := Zeros[j,1];
      end if;
    end for;
    iMax := i;
    if iMax == 0 then
      p := minimumDistance;
      return;
    end if;

    vec := vec[1:iMax];
    vecSorted := Modelica.Math.Vectors.sort(vec);

    // Find p, so that vecSorted[i+1] - vecSorted[i] > 2*minimumDistance
    if vecSorted[1] >= 2*minimumDistance then
       p := minimumDistance;
       return;
    end if;

    i :=1;
    while i <= iMax loop
      if i == iMax then
        p := vecSorted[i] + minimumDistance;
        break;
      elseif vecSorted[i + 1] - vecSorted[i] > 2*minimumDistance then
        p := vecSorted[i] + minimumDistance;
        break;
      else
        i := i + 1;
      end if;
    end while;
  end getReOutsidePolesZeros;

algorithm
  assert(size(B,2)==1 and size(C,1)==1, "System is not a SISO system");
  if size(A,1) == 0 then
     // no states
     gain :=D[1, 1];
     return;
  end if;

  // Determine real value r that is far enough away from all zeros and poles
  r := getReOutsidePolesZeros(Zeros,Poles);

  // With G1(s) = (s-z1)*(s-z2)*...*(s-zm)/((s-p1)*(s-p2)*...(s-pn))"
  // compute kr1 = G1(r)
  (kr1,phi) :=Internal.frequencyEvaluate(1,Zeros,Poles,r,0);
  if abs(phi) > 1 then
     kr1 :=-kr1;
  end if;

  // With G2(s) = C*inv(sI-A)*B + D
  // compute kr2 = G2(r)
  A2 :=-A;
  for i in 1:size(A, 1) loop
     A2[i, i] := A2[i, i] + r;
  end for;
  kr2 :=vector(C)*Modelica.Math.Matrices.solve(A2, vector(B)) + scalar(D);

  // Compute gain
  gain :=kr2/kr1;

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
</html>",
      revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2011-07-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Improved frequency calculation</td>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end frequencyResponseGain;
