within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designInverseDoublePendulumControllerWithObserver
  "Design pole assignment for an inverse double pedulum"

  import Modelica.Utilities.Streams.print;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.StateSpace;
//  input String modelName="Modelica_Controller.Examples.Components.InverseDoublePendulum"  "name of the model to linearize";
  input String modelName="Modelica_LinearSystems2.Controller.Examples.Components.InverseDoublePendulum"
    "name of the model to linearize";

//input Complex pa[6]={-1+0*j, -1+0*j, -6-0.2*j,-6+0.2*j,-6-0.2*j,-6+0.2*j}     "assigned poles";
input Complex pa[6]={Complex(-1,0), Complex(-1,0), Complex(-12,-0.2),Complex(-12,0.2),Complex(-15,-0),Complex(-15,0)}
    "assigned poles";
  input Complex pob[6]=cat(1,fill(-15+0*j,4),fill(-20+0*j,2))
    "assigned observer poles";

 input String fileName=DataDir + "inverseDoublePendulumControllerO.mat"
    "file name for results";

  input Boolean makeAnalysis=false;

 output Real K_pa[:,:] "feedback matrix pole assigment controller";
 output Real M_pa[:,:] "pre filter LQ controller";
 output Real K_ob[:,:] "feedback matrix pole assigment controller";

protected
 input Complex j = Modelica_LinearSystems2.Math.Complex.j();
protected
 Real Q[:,:];
 Real R[:,:];
 Boolean isControllable;
 Boolean isObservable;

  // Determine linear System from Modelica_Controller.Examples.Components.InverseDoublePendulum.mo
 Modelica_LinearSystems2.StateSpace ss = Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);

 Complex p[:]=Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A);
 Modelica_LinearSystems2.StateSpace ss_pa=ss;
// Modelica_LinearSystems2.StateSpace ss_ob=StateSpace(A=transpose(ss.A), B=transpose([1,0,0,0,0,0;0,0,1,0,0,0]), C=transpose(ss.B), D=fill(0,size(ss.B,2),2));
// Modelica_LinearSystems2.StateSpace ssPlant=Modelica_LinearSystems2.StateSpace(A=ss.A, B=ss.B,C=[1,0,0,0,0,0;0,0,1,0,0,0],D=zeros(2,size(ss.B,2)));

// Observer with 3 meassurements
//   Modelica_LinearSystems2.StateSpace ss_ob=StateSpace(A=transpose(ss.A), B=transpose([1,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,1,0]), C=transpose(ss.B), D=fill(0,size(ss.B,2),3));
//   Modelica_LinearSystems2.StateSpace ssPlant=Modelica_LinearSystems2.StateSpace(A=ss.A, B=ss.B,C=[1,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,1,0],D=zeros(3,size(ss.B,2)));

// Observer with 2 meassurements
  Modelica_LinearSystems2.StateSpace ss_ob=StateSpace(A=transpose(ss.A), B=transpose([1,0,0,0,0,0;0,0,1,0,0,0]), C=transpose(ss.B), D=fill(0,size(ss.B,2),2));
  Modelica_LinearSystems2.StateSpace ssPlant=Modelica_LinearSystems2.StateSpace(A=ss.A, B=ss.B,C=[1,0,0,0,0,0;0,0,1,0,0,0],D=zeros(2,size(ss.B,2)), yNames={"s","phi1"}, xNames=ss.xNames, uNames=ss.uNames);

algorithm
 ss.C:=identity(6);
print("The linearized state space system is determined to:\n" + String(ss));
  Modelica_LinearSystems2.Math.Matrices.printMatrix(ss.C,6,"C");
 Modelica_LinearSystems2.Math.Matrices.printMatrix(ss.D,6,"D");

if makeAnalysis then
StateSpace.Analysis.analysis(ssPlant,fileName="inverseDoublePendulum.html");
end if;

isControllable := StateSpace.Analysis.isControllable(ssPlant);

  if isControllable then
      Modelica.Utilities.Streams.print("ssPlant is controllable");
    else
      Modelica.Utilities.Streams.print("ssPlant is not controllable");
    end if;
isObservable := StateSpace.Analysis.isObservable(ssPlant);
  if isObservable then
      Modelica.Utilities.Streams.print("ssPlant is observable");
    else
      Modelica.Utilities.Streams.print("ssPlant is not observable");
    end if;

 print("The eigenvalues are:\n");
Complex.Vectors.print("p",p);

  //####### POLE ASSIGNMENT ##########

  // feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, pa);
  ss_pa.A := ss.A - ss.B*K_pa;

  print("The feedback matrix of the pole assignment controller is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    K_pa,
    6,
    "K_pa"));
  print("eigenvalues of the closed loop system are:\n");
  Modelica_LinearSystems2.Math.Complex.Vectors.print("ev_pa", p);

  writeMatrix(
    fileName,
    "K_pa",
    K_pa,
    true);

  // Pre filter calculation
  M_pa := -Modelica.Math.Matrices.inv([1,0,0,0,0,0]*Matrices.solve2(ss_pa.A,
    ss_pa.B));
  print("Gain for pre filtering:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    M_pa,
    6,
    "M_pa"));
  writeMatrix(
    fileName,
    "M_pa",
    M_pa,
    true);

  K_ob := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss_ob, pob);
  K_ob := transpose(K_ob);

  writeMatrix(
    fileName,
    "nx",
    [size(ssPlant.A,1)],
    true);
 writeMatrix(
    fileName,
    "stateSpace",
    [ssPlant.A,ssPlant.B; ssPlant.C,ssPlant.D],
    true);

  writeMatrix(
    fileName,
    "K_ob2",
    K_ob,
    true);

if makeAnalysis then
StateSpace.Analysis.analysis(ss_pa,fileName="inverseDoublePendulumControlled.html");
end if;
  print("\nok!");
  annotation (interactive=true, Documentation(info="<html>
This example demonstrates how to design pole placement controller to balance an inverted double pendulum. For controller design a linearized model of a (simple) physical system model is used.
The controller is applied to the physical model in Moldelica_Controller library.
<br>
It is also shown how to linearize a modelica model. The linear model is used as a base for control design

</html>"),    Documentation(info="<html>
This example demonstrates how to design a lq-controller or a pole placement controller respectively.
The feedback matrices and a simple pre filter for tracking are save to MATLAB files which can be used in
ModelicaController library.<br>
It is also shown how to linearize a model of a crane trolley system [1]. The linear model is used as a base for control design

<A name=\"References\"><B><FONT SIZE=\"+1\">References</FONT></B></A>
<PRE>
  [1] Föllinger, O. \"Regelungstechnik\", Hüthig-Verlag
</PRE>
</html>"));
end designInverseDoublePendulumControllerWithObserver;
