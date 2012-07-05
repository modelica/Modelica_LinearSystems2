within Modelica_LinearSystems2.Examples;
package DiscreteZerosAndPoles
  "Package with examples to demonstrate the usage of the DiscreteZerosAndPoles record"
  extends Modelica.Icons.ExamplesPackage;
  function poltBode
    "Obsolete function: Bode plot of a DiscreteZerosAndPoles objekt"
    extends Modelica.Icons.Function;
    extends Modelica.Icons.ObsoleteModel;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  protected
    DiscreteZerosAndPoles q=DiscreteZerosAndPoles.q();
    DiscreteZerosAndPoles dzp=0.000160362*(q - 0.818731)*(q + 0.258403)*(q +
        3.59015)*(q^2 - 1.80063*q + 0.818731)/((q - 0.904837)*(q - 0.860708)*(q
        ^2 - 1.89533*q + 0.904837)*(q^2 - 1.79161*q + 0.818731));
    Real Ts=0.1;
    Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
  algorithm
    dzp.Ts := Ts;
    dzp.method := method;
    DiscreteZerosAndPoles.Plot.bode(dzp);

  annotation (interactive=true, Documentation(info="<html>
<p>This function is obsolete. Use <a href=\"Modelica://Modelica_LinearSystems2.Examples.DiscreteZerosAndPoles.plotBode\">plotBode</a> instead.</p>
</html>"));
  end poltBode;

  function plotBode "Bode plot of a DiscreteZerosAndPoles objekt"
    extends Modelica.Icons.Function;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

  protected
    DiscreteZerosAndPoles q=DiscreteZerosAndPoles.q();
    DiscreteZerosAndPoles dzp=0.000160362*(q - 0.818731)*(q + 0.258403)*(q +
        3.59015)*(q^2 - 1.80063*q + 0.818731)/((q - 0.904837)*(q - 0.860708)*(q
        ^2 - 1.89533*q + 0.904837)*(q^2 - 1.79161*q + 0.818731));
    Real Ts=0.1;
    Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.StepExact;
  algorithm
    dzp.Ts := Ts;
    dzp.method := method;
    DiscreteZerosAndPoles.Plot.bode(dzp);

  annotation (interactive=true);
  end plotBode;
end DiscreteZerosAndPoles;
