within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotStep "Step plot example"

  import Modelica_LinearSystems2.ZerosAndPoles;

  annotation (interactive=true, Documentation(info="<html>
<p>
Computes the impulse response of the system tf =1/s^2 + s + 1.
</html>"));

protected
  ZerosAndPoles zp=ZerosAndPoles(k=1,n1=fill(0,0),d2=[1,1]);

algorithm
 Modelica_LinearSystems2.ZerosAndPoles.Plot.step(        zp=zp, dt=0.1, tSpan=10);

end plotStep;
