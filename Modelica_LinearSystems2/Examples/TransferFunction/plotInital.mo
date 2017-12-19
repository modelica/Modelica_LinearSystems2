within Modelica_LinearSystems2.Examples.TransferFunction;
function plotInital "Example plotting initial condition response"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.TransferFunction;

  input TransferFunction tf=TransferFunction({1}, {1,1,1});

protected
  Real y0=1 "Initial state vector";
algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.initialResponse(tf=tf, y0=y0, dt=0.1, tSpan=10);

  annotation (
    __Dymola_interactive=true,
    Documentation(info="<html>
<p>
This example computes the initial condition response of the of the system
<i>tf&nbsp;=&nbsp;1/s^2&nbsp;+&nbsp;s&nbsp;+&nbsp;1</i>
to the initial condition <i>y0=y0, dt=0.1, tSpan=10</i>.
</p>
</html>
"));
end plotInital;
