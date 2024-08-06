within Modelica_LinearSystems2.Examples.StateSpace;
function plotTimeResponse "Time response plot example"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.Utilities.Types.TimeResponse response=Modelica_LinearSystems2.Utilities.Types.TimeResponse.Step "Type of time response";

  input Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
    A=[-1,1; 0,-2],
    B=[4,0; 0,2],
    C=[2,0; 0,3; -1, 2],
    D=[0,0; 0,0; 1, 0],
    yNames={"phi1", "w", "a"},
    uNames={"tau1", "f1"});

algorithm
  Modelica_LinearSystems2.StateSpace.Plot.timeResponse(ss=ss);

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
Computes the impulse response of the system
StateSpace <em>sc = StateSpace(A=[-1,1;0,-2],B=[1, 0;0, 1],C=[1,0; 0,1],D=[0, 0; 0, 0])</em>.
</p>
</html>"));
end plotTimeResponse;
