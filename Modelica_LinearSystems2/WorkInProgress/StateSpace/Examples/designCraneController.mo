within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designCraneController
  "Design pole assignment and LQ controller for an overhead crane"
  import Modelica.Utilities.Streams;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.ComplexMathAdds;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input String modelName=
      "Modelica_Controller.Examples.Components.Pendulum_small"
    "name of the model to linearize";

  input Complex pa[4]={-3.5 + 0*j,-3.5 + 0*j,-3.5 - 0.5*j,-3.5 + 0.5*j}
    "assigned poles";

  output Real K_lq[:,:] "feedback matrix LQ controller";
  output Real K_pa[:,:] "feedback matrix pole assignment controller";
  output Real M_lq[:,:] "pre filter LQ controller";
  output Real M_pa[:,:] "pre filter pole assignment controller";
protected
  input Complex j = Modelica.ComplexMath.j;

protected
  Real Q[:,:];
  Real R[:,:];

// Determine linear System from Modelica_Controller.Examples.Pendulum.mo
  Modelica_LinearSystems2.StateSpace ss=
      Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);
//Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(A= [0,1,0,0;0,0,39.24,0;0,0,0,1;0,0,-4.904,0], B=[0;1e-3;0;-1e-4], C=[1,0,1,0], D=[0]);

  Complex p[:]=ComplexMathAdds.eigenValues(ss.A);
  Modelica_LinearSystems2.StateSpace ss_lq=ss;
  Modelica_LinearSystems2.StateSpace ss_pa=ss;

algorithm
  Streams.print("The linearized state space system is determined to:\n" + String(ss));

// ####### LQ CONTROLLER #######

// eigenvalues of open loop system
  p := ComplexMathAdds.eigenValues(ss.A);
  Streams.print("eigenvalues of the open loop system are:\n");
  ComplexMathAdds.Vectors.print("ev", p);

// Calculate feesback matrix of a lq controller (Riccati) with weighting matrices Q and R
  Q := [1,0,300,1000; 0,100,0,3000; 0,0,1000,0; 0,10,0,1];
  Q := Q*transpose(Q);
  R := identity(size(ss.B, 2));
  (K_lq) := Modelica_LinearSystems2.StateSpace.Design.lqr(
    ss,
    Q,
    R);

  Streams.print("\nXXXXXXXXXXX:\n");

  ss_lq.A := ss.A - ss.B*K_lq;
  Streams.print("The feedback matrix of the lq controller is:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    K_lq,
    6,
    "K_lq"));
  Streams.writeRealMatrix(
    DataDir + "craneController_small.mat",
    "K_lq",
    K_lq);

// eigenvalues of closed loop system
  p := ComplexMathAdds.eigenValues(ss_lq.A);
  Streams.print("eigenvalues of the closed loop system are:\n");
  ComplexMathAdds.Vectors.print("ev_lq", p);

// Pre filter calculation
  M_lq := -MatricesMSL.inv([1,0,0,0]*MatricesMSL.solve2(ss_lq.A, ss_lq.B));
  Streams.print("Gain for pre filtering:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    M_lq,
    6,
    "M_lq"));
  Streams.writeRealMatrix(
    DataDir + "craneController_small.mat",
    "M_lq",
    M_lq,
    true);

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
  ComplexMathAdds.Vectors.print("ev_pa", p);

  Streams.writeRealMatrix(
    DataDir + "craneController_small.mat",
    "K_pa",
    K_pa,
    true);

// Pre filter calculation
  M_pa := -MatricesMSL.inv([1,0,0,0]*MatricesMSL.solve2(ss_pa.A, ss_pa.B));
  Streams.print("Gain for pre filtering:\n" +
    Modelica_LinearSystems2.Math.Matrices.printMatrix(
    M_pa,
    6,
    "M_pa"));
  Streams.writeRealMatrix(
    DataDir + "craneController_small.mat",
    "M_pa",
    M_pa,
    true);

  Streams.print("\nok!");
annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates how to design a lq-controller or a pole placement controller
respectively. Compared with example <strong>craneController</strong>, the plant is smaller to
achieve suitable dynamics for animation.
The feedback matrices and a simple pre filter for tracking are save to MATLAB files
which can be used in ModelicaController library.
</p>
<p>
It is also shown how to linearize a model of a crane trolley system [1].
The linear model is used as a base for control design
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] F&ouml;llinger O.:</dt>
<dd> <strong>Regelungstechnik</strong>.
     H&uuml;thig-Verlag.<br>&nbsp;</dd>
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
end designCraneController;
