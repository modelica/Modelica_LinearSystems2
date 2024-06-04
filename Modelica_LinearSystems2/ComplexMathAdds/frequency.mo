within Modelica_LinearSystems2.ComplexMathAdds;
function frequency "Frequency and damping of conjugated complex pole pair"
  extends Modelica.Icons.Function;

  input Complex c "Complex number";
  output Modelica.Units.SI.Frequency f "Frequency of c (= c.im in Hz)";
  output Real damping "Damping of c (= c.re/c.im)";

protected
  Real abs_ev = sqrt(c.re^2 + c.im^2);
algorithm
  f := if abs(c.im) > 10*Modelica.Constants.eps then abs_ev/(2*Modelica.Constants.pi) else 0;
  damping := if abs(c.im) > 10*Modelica.Constants.eps then
    if abs_ev > Modelica.Constants.eps then -c.re/abs_ev else 0.0 else 1.0;
  annotation(Inline=true);
end frequency;
