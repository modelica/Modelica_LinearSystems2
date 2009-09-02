within Modelica_LinearSystems2.Examples.StateSpace;
function transformationExtract
  "Example how to extract input/output related subsystems from state space system record"

protected
  Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
      A=[-1,1; 0,-2],
      B=[1,0; 0,1],
      C=[1,0; 0,1],
      D=[0,0; 0,0]);
Integer i1[1]={1};
Integer i2[1]={1};

Modelica_LinearSystems2.StateSpace subSys(
    redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
    redeclare Real B[size(ss.B, 1),size(i1,1)],
    redeclare Real C[size(i2,1),size(ss.C, 2)],
    redeclare Real D[size(i2, 1),size(i1, 1)]) "Subsystem state space record";
algorithm
subSys :=  Modelica_LinearSystems2.StateSpace.Transformation.extract(
    ss,
    i1,
    i2);
 Modelica.Utilities.Streams.print("Subsystem is " +String(subSys));
end transformationExtract;
