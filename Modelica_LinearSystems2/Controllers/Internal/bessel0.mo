within Modelica_LinearSystems2.Controllers.Internal;
function bessel0
  "Polynomial approximation of the zeroth order modified Bessel function"
  extends Modelica.Icons.Function;

  input Real x;
  output Real y;
protected
  Real ax;
  Real a;

algorithm
  ax := abs(x);
  if ax < 3.75 then
    a := (x/3.75)^2;
    y := 1 + a*(3.5156229 + a*(3.0899424 + a*(1.2067492 + a*(0.2659732 + a*
      (0.0360768 + a*0.0045813)))));
  else
    a := 3.75/ax;
    y := exp(ax)/sqrt(ax)*(0.39894228 + a*(0.01328592 + a*(0.00225319 + a*(
      -0.00157565 + a*(0.00916281 + a*(-0.02057706 + a*(0.02635537 + a*(-0.01647633
       + a*0.00392377))))))));
  end if;
  annotation (
    Documentation(info="<html>
<p>
Polynomial approximation of the zeroth order modified Bessel function.
The algorithm is taken from&nbsp;[1].
The function is used to calculate the Kaiser-window via <em>calcWindow</em>.
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>[1] H. W. Press, S.A. Teukolsky, W. Vetterling:
<dd><strong>Numerical Reciepes in C: The Art of Scientific Computing</strong><br>
       Cambridge UP, 1988
</dl>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2002-07-10</td>
    <td valign=\"top\">Nico Walther</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end bessel0;
