within Modelica_LinearSystems2.Controller.Types;
type PID_representation = enumeration(
    gains "gains representation: G_PID = kp + ki/s + kd*s",
    timeConstants
      "time constants representation: G_PID = k*(1 + 1/Ti/s + Td*s)")
  "Enumeration defining the representation of a PID controller" 
    annotation (Evaluate=true, Documentation(info="<html>
 
</html>"));
