within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designCraneControllerWithObserver
  "Design pole assignment controller and observer for an overhead crane"
  import Modelica.Utilities.Streams.print;
  import Modelica.Utilities.Streams.writeRealMatrix;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Complex;

  input String modelName="Modelica_Controller.Examples.Components.Pendulum_small"
    "Name of the model to linearize";
  input Complex pa[4]={-3.5+0*j,-3.5+0*j, -3.5-0.5*j, -3.5+0.5*j}
    "Assigned system poles";
  input Complex pob[4]={-10+0*j,-10+0*j,-10+0*j,-10+0*j}
    "Assigned observer poles";
  input String fileName=DataDir + "craneController_small.mat"
    "File name for results";

protected
  input Complex j = Modelica.ComplexMath.j;
public
  output Real K_ob[:,:] "Feedback matrix pole assignment controller";
  output Real K_pa[:,:] "Feedback matrix pole assignment controller";
  output Real M_pa[:,:] "Pre filter LQ controller";

// Determine linear System from Modelica_Controller.Examples.Pendulum.mo
protected
  Modelica_LinearSystems2.StateSpace ss=
  Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);

  Complex p[:];//=Modelica_LinearSystems2.Math.Complex.eigenValues(ss.A);
  Modelica_LinearSystems2.StateSpace ss_pa=ss;
  Modelica_LinearSystems2.StateSpace ss_ob=StateSpace(A=transpose(ss.A), B=[1,0;0,1;0,0;0,0], C=transpose(ss.B), D=[0,0]);
  Modelica_LinearSystems2.StateSpace ssPlant=StateSpace(A=ss.A, B=ss.B,C=[1,0,0,0;0,1,0,0], D=[0;0]);

algorithm
  print("The linearized state space system is determined to:\n" + String(ssPlant));

 // eigenvalues of open loop system
  p :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(ss.A);
  print("eigenvalues of the open loop system are:\n");
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("ev", p);

//####### POLE ASSIGNMENT ##########

// feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, pa);
  ss_pa.A := ss.A - ss.B*K_pa;
  print("The feedback matrix of the pole assignment controller is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(K_pa, 6, "K_pa"));
  print("eigenvalues of the closed loop system are:\n");
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("ev_pa", p);

  writeRealMatrix(fileName, "K_pa", K_pa, false);

// Pre filter calculation
  M_pa := -Modelica.Math.Matrices.inv([1,0,0,0]*Matrices.solve2(ss_pa.A, ss_pa.B));
  print("Gain for pre filtering:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(M_pa, 6, "M_pa"));
  writeRealMatrix(fileName, "M_pa", M_pa, true);

// observer feedback
  (K_ob,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss_ob, pob);
  K_ob := transpose(K_ob);

  writeRealMatrix(fileName, "stateSpace", [ssPlant.A, ssPlant.B; ssPlant.C, ssPlant.D], true);
// write matrix dimension nx
  writeRealMatrix(fileName, "nx", [size(ssPlant.A, 1)], true);
  print("The feedback matrix of the observer system is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(K_ob, 6, "K_ob"));
  ss_ob.A := ss.A - K_ob*ssPlant.C;

  print("eigenvalues of the observer system are:\n");
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("ev_pob", pob);
  writeRealMatrix(fileName, "K_ob", K_ob, true);

  print("\nok!");
  annotation (
    Documentation(info="<html>
<p>
This example demonstrates how to use pole placement function assignPolesMI to
design a pole placement controller and an observer also using assignPolesMI as well.
Function assignPoles is based on the Schur form of the state space matrix an also
allows partial poles shifting [1].
</p>
<p>
An observer can be designed by applying assignPolesMI to the dual system (A', C', B').
The results, which are controller feedback matrix, observer feedback matrix, system
model for usage in observer and a simple pre filter for tracking are saved to MATLAB
files which can be used in ModelicaController library.
It is also shown how to linearize a model of a crane trolley system [1].
The linear model is then used as a base for controller design.
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Varga A. (1981):</dt>
<dd> <b>A Schur method for pole assignment</b>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>
</html>"));
end designCraneControllerWithObserver;
