within Modelica_LinearSystems2.Examples.TransferFunction;
function plotImpulse "Example plotting impulse response"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  input TransferFunction tf=TransferFunction(n={1}, d={1,1,1});

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.impulse(
    tf=tf,
    dt=0.1,
    tSpan=10);

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This example computes the impulse response of the system
<i>tf&nbsp;=&nbsp;1/s^2&nbsp;+&nbsp;s&nbsp;+&nbsp;1</i>.
</p>
</html>"));
end plotImpulse;
