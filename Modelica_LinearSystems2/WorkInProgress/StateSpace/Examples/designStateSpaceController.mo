within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designStateSpaceController
  "Demonstration of controller design for a state space system"
  import Modelica.Utilities.Streams.print;
  import Modelica.ComplexMath;
  import Modelica.Math.Matrices;
  import Modelica_LinearSystems2.ComplexMathAdds;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.StateSpace;

  input StateSpace ss=StateSpace(
      A=[0,1,0,0; 0,0,39.2,0; 0,0,0,1; 0,0,-49.0,0],
      B=[1,1,0,0;0,1,1,0;0,0,1,1;1,1,1,1],
      C=identity(4),
      D=zeros(4,4),
      yNames={"s","v","phi","w"},
      uNames={"f","u2","u3","u4"},
      xNames={"s","v","phi","w"});

   input Complex pa[:]={Complex(-1,0),Complex(-2,0),Complex(-2, -0.2),Complex(-2,0.2)}
    "assigned poles";

  output Real K_pa[:,:] "feedback matrix pole assignment controller";
  output Real M_pa[:,:] "pre filter pole assignment controller";
  output Complex po[size(ss.A,1)];

protected
  Complex assignedPoles[:]=pa;
  StateSpace ss_pa=ss;

algorithm
//####### POLE ASSIGNMENT ##########
// return the Schur form representation
  (,,po) := StateSpace.Design.assignPolesMI2(ss,assignedPoles,size(assignedPoles,1),true);
  ComplexMathAdds.Vectors.print("poSchur", po);
  StateSpace.Analysis.analysis(ss,fileName="openloopSystem.html");

// feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,po) := StateSpace.Design.assignPolesMI2(ss, assignedPoles);
  ss_pa.A := ss.A - ss.B*K_pa;
  StateSpace.Analysis.analysis(ss_pa,fileName="closedloopSystem.html");
  print("The feedback matrix of the pole assignment controller is:\n" +
    Matrices.toString(K_pa, "K_pa", 6));
  print("eigenvalues of the closed loop system are:\n");
  ComplexMathAdds.Vectors.print("ev_pa", po);
// Pre filter calculation
//   M_pa := -Modelica.Math.Matrices.inv([1,0,0,0;0,0,1,0]*Modelica_LinearSystems2.Math.Matrices.solve2(ss_pa.A, ss_pa.B));
//   print("Gain for pre filtering:\n" +
//     Matrices.toString(M_pa, "M_pa", 6));

  assignedPoles := ComplexMath.Vectors.reverse(assignedPoles);
  assignedPoles[3]:=Complex(-1);
  assignedPoles[4]:=Complex(-2);
// feedback matrix of a pole assignment controller with assigned eigenvalues pa
  (K_pa,,po) := StateSpace.Design.assignPolesMI2(ss, assignedPoles);
  ss_pa.A := ss.A - ss.B*K_pa;
  StateSpace.Analysis.analysis(ss_pa,fileName="closedloopSystem2.html");
  print("The feedback matrix of the pole assignment controller is:\n" +
    Matrices.toString(K_pa, "K_pa", 6));
  print("eigenvalues of the closed loop system are:\n");
  ComplexMathAdds.Vectors.print("ev_pa", po);
// Pre filter calculation
//   M_pa := -Modelica.Math.Matrices.inv([1,0,0,0;0,0,1,0]*Modelica_LinearSystems2.Math.Matrices.solve2(ss_pa.A, ss_pa.B));
//   print("Gain for pre filtering:\n" +
//     Matrices.toString(M_pa, "M_pa", 6));

  print("\nok!");
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates how to design a lq-controller or a pole placement controller respectively.
Compared with example <strong>craneController</strong>,
the plant is smaller to achieve suitable dynamics for animation.
The feedback matrices and a simple pre filter for tracking are save to MATLAB files which can be used in
ModelicaController library.
</p>
<p>
It is also shown how to linearize a model of a crane trolley system [1]. The linear model is used as a base for control design
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] F&ouml;llinger O.:</dt>
<dd> <strong>Regelungstechnik</strong>.
     H&uuml;thig-Verlag.<br>&nbsp;</dd>
</dl>

</html>"));
end designStateSpaceController;
