within Modelica_LinearSystems2.Examples.DiscreteZerosAndPoles;
function plotBode
  "Example for construction of DiscreteZerosAndPoles objekt and plot of the Bode diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.DiscreteZerosAndPoles;

protected
  DiscreteZerosAndPoles q=DiscreteZerosAndPoles.q();
  DiscreteZerosAndPoles dzp=0.000160362*(q - 0.818731)*(q + 0.258403)*(q +
    3.59015)*(q^2 - 1.80063*q + 0.818731)/((q - 0.904837)*(q - 0.860708)*(q
    ^2 - 1.89533*q + 0.904837)*(q^2 - 1.79161*q + 0.818731));
  Real Ts=0.1;
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact;
algorithm
  dzp.Ts := Ts;
  dzp.method := method;
  DiscreteZerosAndPoles.Plot.bode(dzp);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example shows how to construct a discrete zeros and poles system and to plot the Bode diagram
with automatic determination of the frequency range to plot.
</p>
</html>
"));
end plotBode;
