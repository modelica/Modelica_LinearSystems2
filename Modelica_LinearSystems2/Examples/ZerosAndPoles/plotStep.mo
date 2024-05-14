within Modelica_LinearSystems2.Examples.ZerosAndPoles;
function plotStep "Step plot example"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.ZerosAndPoles;

protected
  ZerosAndPoles zp = ZerosAndPoles(k=1,n1=fill(0,0),d2=[1,1]);

algorithm
  ZerosAndPoles.Plot.step(zp=zp, dt=0.1, tSpan=10);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the impulse response of the system tf =1/s^2 + s + 1.
</p>
</html>"));
end plotStep;
