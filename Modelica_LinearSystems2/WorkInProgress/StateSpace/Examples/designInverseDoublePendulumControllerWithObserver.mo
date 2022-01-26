within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designInverseDoublePendulumControllerWithObserver
  "Design pole assignment for an inverse double pedulum"

  import Modelica.Utilities.Streams;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica.ComplexMath.j;
  import Modelica_LinearSystems2.ComplexMathAdds;
  import Modelica_LinearSystems2;
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

  output Real K_pa[:,:] "feedback matrix pole assignment controller";
  output Real M_pa[:,:] "pre filter LQ controller";
  output Real K_ob[:,:] "feedback matrix pole assignment controller";

protected
 Real Q[:,:];
 Real R[:,:];
 Boolean isControllable;
 Boolean isObservable;

  // Determine linear System from Modelica_Controller.Examples.Components.InverseDoublePendulum.mo
 Modelica_LinearSystems2.StateSpace ss = Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);

 Complex p[:]=ComplexMathAdds.eigenValues(ss.A);
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
  Streams.print("The linearized state space system is determined to:\n" + String(ss));
  Matrices.printMatrix(ss.C,6,"C");
  Matrices.printMatrix(ss.D,6,"D");

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

  Streams.print("The eigenvalues are:\n");
  Modelica_LinearSystems2.ComplexMathAdds.Vectors.print("p",p);

  //####### POLE ASSIGNMENT ##########

  // feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, pa);
  ss_pa.A := ss.A - ss.B*K_pa;

  Streams.print("The feedback matrix of the pole assignment controller is:\n" +
    Matrices.printMatrix(
    K_pa,
    6,
    "K_pa"));
  Streams.print("eigenvalues of the closed loop system are:\n");
  ComplexMathAdds.Vectors.print("ev_pa", p);

  Streams.writeRealMatrix(
    fileName,
    "K_pa",
    K_pa,
    true);

  // Pre filter calculation
  M_pa := -MatricesMSL.inv([1,0,0,0,0,0]*MatricesMSL.solve2(ss_pa.A,
    ss_pa.B));
  Streams.print("Gain for pre filtering:\n" +
    Matrices.printMatrix(
    M_pa,
    6,
    "M_pa"));
  Streams.writeRealMatrix(
    fileName,
    "M_pa",
    M_pa,
    true);

  K_ob := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss_ob, pob);
  K_ob := transpose(K_ob);

  Streams.writeRealMatrix(
    fileName,
    "nx",
    [size(ssPlant.A,1)],
    true);
  Streams.writeRealMatrix(
    fileName,
    "stateSpace",
    [ssPlant.A,ssPlant.B; ssPlant.C,ssPlant.D],
    true);

  Streams.writeRealMatrix(
    fileName,
    "K_ob2",
    K_ob,
    true);

  if makeAnalysis then
    StateSpace.Analysis.analysis(ss_pa,fileName="inverseDoublePendulumControlled.html");
  end if;
  Streams.print("\nok!");
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates how to design pole placement controller to balance an inverted double pendulum. For controller design a linearized model of a (simple) physical system model is used.
The controller is applied to the physical model in Moldelica_Controller library.
</p>
<p>
It is also shown how to linearize a modelica model. The linear model is used as a base for control design
</p>
</html>"));
end designInverseDoublePendulumControllerWithObserver;
