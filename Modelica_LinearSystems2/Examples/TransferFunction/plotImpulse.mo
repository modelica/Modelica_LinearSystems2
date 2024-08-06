within Modelica_LinearSystems2.Examples.TransferFunction;
function plotImpulse "Impulse plot example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  input TransferFunction tf=TransferFunction(n={1}, d={1,1,1});

algorithm
  TransferFunction.Plot.impulse(tf, 0.1, 10);

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
Computes the impulse response of the system tf =1/s^2 + s + 1.
</p>
</html>"));
end plotImpulse;
