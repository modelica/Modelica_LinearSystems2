within Modelica_LinearSystems2.Examples.TransferFunction;
function plotStep "Step plot example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  input TransferFunction tf=TransferFunction({1}, {1,1,1});

algorithm
  TransferFunction.Plot.step(tf, 0.1, 10);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the step response of the system tf =1/s^2 + s + 1.
</p>
</html>"));
end plotStep;
