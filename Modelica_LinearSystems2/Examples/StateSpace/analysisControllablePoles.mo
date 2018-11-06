within Modelica_LinearSystems2.Examples.StateSpace;
function analysisControllablePoles
  "Example to check controllability of a state space system and print the controllable poles"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica.Utilities.Streams.print;

  input StateSpace ssi=Modelica_LinearSystems2.StateSpace(
    A=[1,0,0,0,0,0; 1,0,0,0,0,0; 0,2,3,0,78,6; 1,1,2,2,3,3; 10,13,34,0,0,1; 0,
       0,0,2,0,0],
    B=[0,0; 0,0; 0,0; 0,0; 1,0; 0,0],
    C=[1,0,1,0,1,0; 0,1,0,1,0,1; 0,1,0,1,0,1],
    D=[0,0; 0,0; 0,0]);

  input Boolean systemOnFile=false
    "True, if state space system is defined on file"
    annotation(Dialog(group="system data definition"),choices(checkBox=true));
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored"
    annotation(Dialog(group="system data definition",loadSelector(filter="MAT files (*.mat);; All files (*.*)",
     caption="state space system data file"),enable = systemOnFile));
  input String matrixName="ABCD" "Name of the state space system matrix"
    annotation(Dialog(group="system data definition",enable = systemOnFile));

  output Boolean ok;

protected
  StateSpace ss = if systemOnFile then
    Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else ssi;
  Real cPoles[:,2] "controllable poles";
  Real ncPoles[:,2] "uncontrollable poles";

algorithm
  (cPoles,ncPoles) := StateSpace.Internal.controllablePoles(ss);
  if size(ncPoles, 1) == 0 then
    print("\nThe system is controllable\nThe poles are" +
      Modelica_LinearSystems2.Math.Matrices.printMatrix(cPoles,6,""));
  else
    print("\nThe system is not controllable\nThe uncontrollable poles are" +
      Modelica_LinearSystems2.Math.Matrices.printMatrix(ncPoles,6,""));
  end if;

  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example shows the usage of function 
<a href=\"modelica://Modelica_LinearSystems2.StateSpace.Internal.controllablePoles\">StateSpace.Internal.controllablePoles</a>
which is to compute the controllable and uncontrollable poles of a state space system.
</p>
</html>"));
end analysisControllablePoles;
