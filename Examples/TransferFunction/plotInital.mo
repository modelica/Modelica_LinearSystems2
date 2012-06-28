within Modelica_LinearSystems2.Examples.TransferFunction;
function plotInital "Initial condition plot example"

  import Modelica_LinearSystems2.TransferFunction;

 input TransferFunction tf=TransferFunction({1}, {1,1,1});

protected
  Real y0=1 "Initial state vector";
algorithm
  Modelica_LinearSystems2.TransferFunction.Plot.initialResponse(tf=tf, y0=y0, dt=0.1, tSpan=10);

  annotation (interactive=true, Documentation(info="<html>

</html>"));
end plotInital;
