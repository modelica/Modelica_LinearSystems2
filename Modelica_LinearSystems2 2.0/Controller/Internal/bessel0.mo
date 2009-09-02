within Modelica_LinearSystems2.Controller.Internal;
function bessel0
  "Polynomial approximation of the zeroth order modified Bessel function"

  input Real x;
  output Real y;
protected
  Real ax;
  Real a;
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Documentation(info="<HTML>
<p>
Polynomial approximation of the zeroth order modified Bessel function.
The algorithm is taken from
</p>
<dl>
<dt>H. W. Press, S.A. Teukolsky, W. Vetterling:
<dd><b>Numerical Reciepes in C: The Art of Scientific Computing</b><br>
       Cambridge UP, 1988
</dl>
<p>
The function is used to calculate the Kaiser-window via
<i>calcWindow</i>.
</p>
<p><b>Release Notes:</b></p>
<ul>
<li><i>July 10, 2002</i>
       by Nico Walther<br>
       Realized.</li>
</ul>
</HTML>
"), Window(
      x=0.44,
      y=0.25,
      width=0.49,
      height=0.49));
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
end bessel0;
