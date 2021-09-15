within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designCraneControllerWithObserver
  "Design pole assignment controller and observer for an overhead crane"
  import Modelica.Utilities.Streams;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;

  input String modelName="Modelica_Controller.Examples.Components.Pendulum_small"
    "name of the model to linearize";

  input Complex pa[4]={-3.5+0*j,-3.5+0*j, -3.5-0.5*j, -3.5+0.5*j}
    "assigned system poles";

  input Complex pob[4]={-10+0*j,-10+0*j,-10+0*j,-10+0*j}
    "assigned observer poles";
 input String fileName=DataDir + "craneController_small.mat"
    "file name for results";

protected
  input Complex j = Modelica.ComplexMath.j;
public
  output Real K_ob[:,:] "feedback matrix pole assignment controller";
  output Real K_pa[:,:] "feedback matrix pole assignment controller";
  output Real M_pa[:,:] "pre filter LQ controller";

// Determine linear System from Modelica_Controller.Examples.Pendulum.mo
protected
  Modelica_LinearSystems2.StateSpace ss=
  Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);

  Complex p[:];//=Modelica_LinearSystems2.ComplexMathAdds.eigenValues(ss.A);
  Modelica_LinearSystems2.StateSpace ss_pa=ss;
  Modelica_LinearSystems2.StateSpace ss_ob=StateSpace(A=transpose(ss.A), B=[1,0;0,1;0,0;0,0], C=transpose(ss.B), D=[0,0]);
  Modelica_LinearSystems2.StateSpace ssPlant=StateSpace(A=ss.A, B=ss.B,C=[1,0,0,0;0,1,0,0], D=[0;0]);

algorithm
  Streams.print("The linearized state space system is determined to:\n" + String(ssPlant));

 // eigenvalues of open loop system
  p := Modelica_LinearSystems2.ComplexMathAdds.eigenValues(ss.A);
  Streams.print("eigenvalues of the open loop system are:\n");
  Modelica_LinearSystems2.ComplexMathAdds.Vectors.print("ev", p);

//####### POLE ASSIGNMENT ##########

// feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss, pa);
  ss_pa.A := ss.A - ss.B*K_pa;
  Streams.print("The feedback matrix of the pole assignment controller is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    K_pa,
    6,
    "K_pa"));
  Streams.print("eigenvalues of the closed loop system are:\n");
  Modelica_LinearSystems2.ComplexMathAdds.Vectors.print("ev_pa", p);

  Streams.writeRealMatrix(
    fileName,
    "K_pa",
    K_pa,
    false);

// Pre filter calculation
  M_pa := -Modelica.Math.Matrices.inv([1,0,0,0]*Modelica.Math.Matrices.solve2(ss_pa.A, ss_pa.B));
  Streams.print("Gain for pre filtering:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    M_pa,
    6,
    "M_pa"));
  Streams.writeRealMatrix(
    fileName,
    "M_pa",
    M_pa,
    true);

// observer feedback
  (K_ob,,p) := Modelica_LinearSystems2.StateSpace.Design.assignPolesMI(ss_ob, pob);
  K_ob := transpose(K_ob);

  Streams.writeRealMatrix(
    fileName,
    "stateSpace",
    [ssPlant.A,ssPlant.B; ssPlant.C,ssPlant.D],
    true);
// write matrix dimension nx
  Streams.writeRealMatrix(
    fileName,
    "nx",
    [size(ssPlant.A,1)],
    true);
  Streams.print("The feedback matrix of the observer system is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    K_ob,
    6,
    "K_ob"));
  ss_ob.A := ss.A - K_ob*ssPlant.C;

  Streams.print("eigenvalues of the observer system are:\n");
  Modelica_LinearSystems2.ComplexMathAdds.Vectors.print("ev_pob", pob);
  Streams.writeRealMatrix(
    fileName,
    "K_ob",
    K_ob,
    true);

  Streams.print("\nok!");
  annotation (__Dymola_interactive=true, Documentation(info="<html>
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
<dd> <strong>A Schur method for pole assignment</strong>.
     IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519.<br>&nbsp;</dd>
</dl>
</html>"),    Documentation(info="<html>
This example demonstrates how to design a lq-controller or a pole placement controller respectively.
The feedback matrices and a simple pre filter for tracking are save to MATLAB files which can be used in
ModelicaController library.<br>
It is also shown how to linearize a model of a crane trolley system [1]. The linear model is used as a base for control design

<A name=\"References\"><B><FONT SIZE=\"+1\">References</FONT></B></A>
<PRE>
  [1] F&ouml;llinger, O. &quot;Regelungstechnik&quot;, H&uuml;thig-Verlag
</PRE>
</html>"));
end designCraneControllerWithObserver;
