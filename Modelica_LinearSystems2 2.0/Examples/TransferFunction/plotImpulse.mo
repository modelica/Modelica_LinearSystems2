within Modelica_LinearSystems2.Examples.TransferFunction;
function plotImpulse "Impulse plot example"

  import Modelica_LinearSystems2.TransferFunction;

  annotation (interactive=true, Documentation(info="<html>
<p>
Computes the impulse response of the system tf =1/s^2 + s + 1.
</html>"));

 input TransferFunction tf=TransferFunction(n={1}, d={1,1,1});

algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.impulse(
    tf=tf,
    dt=0.1,
    tSpan=10);

end plotImpulse;
