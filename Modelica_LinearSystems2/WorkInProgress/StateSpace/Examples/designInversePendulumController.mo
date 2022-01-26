within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function designInversePendulumController
  "Design pole assignment for an inverse pedulum"

  import Modelica.Utilities.Streams;
  import MatricesMSL = Modelica.Math.Matrices;
  import Modelica.ComplexMath.j;
  import Modelica_LinearSystems2.ComplexMathAdds;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices;

  input String modelName="Modelica_Controller.Examples.Components.InversePendulum_small"
    "Name of the model to linearize";
  input Complex pa[4]={-5+0*j,-5+0*j,-5.0-0.25*j,-5.0+0.25*j} "Assigned poles";

  input String fileName=DataDir + "inversePendulumController_small.mat"
    "File name for results";

public
  output Real K_pa[:,:] "feedback matrix pole assignment controller";
  output Real M_pa[:,:] "pre filter LQ controller";
// Determine linear System from Modelica_Controller.Examples.Pendulum.mo
protected
  Modelica_LinearSystems2.StateSpace ss = Modelica_LinearSystems2.StateSpace.Import.fromModel(modelName);

  Complex p[:]=ComplexMathAdds.eigenValues(ss.A);

  Modelica_LinearSystems2.StateSpace ss_pa=ss;

algorithm
  Streams.print("The linearized state space system is determined to:\n" +String(ss));

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
  M_pa := -MatricesMSL.inv([1,0,0,0]*MatricesMSL.solve2(ss_pa.A, ss_pa.B));
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

  Streams.print("\nok!");

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This example demonstrates how to design pole placement controller to balance an inverted pendulum. For controller design a linearized model of a (simple) physical system model is used.
The controller is applied to the physical model in Moldelica_Controller library.
</p>
<p>
It is also shown how to linearize a modelica model. The linear model is used as a base for control design
</p>
</html>"));
end designInversePendulumController;
