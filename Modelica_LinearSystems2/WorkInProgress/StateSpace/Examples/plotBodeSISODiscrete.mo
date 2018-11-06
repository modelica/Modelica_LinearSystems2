within Modelica_LinearSystems2.WorkInProgress.StateSpace.Examples;
function plotBodeSISODiscrete
  "Plots the Bode diagram of continuous state space system and the corresponding discrete state space system and  with automatic determination of the frequency range to plot "

  import Modelica;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace;
  import Modelica_LinearSystems2.TransferFunction;

  input Boolean systemOnFile=false
    "True, if state space system is defined on file";
  input String fileName="NoName" "file where matrix [A, B; C, D] is stored";

  input Real A[:,size(A, 1)]=[-1.0,0.0,0.0; 0.0,-2.0,0.0; 0.0,0.0,-3.0];
  input Real B[size(A, 2),:]=[0.0,1.0; 1.0,1.0; -1.0,0.0];
  input Real C[:,size(A, 1)]=[0.0,1.0,1.0; 1.0,1.0,1.0];
  input Real D[size(C, 1),size(B, 2)]=[1.0,0.0; 0.0,1.0];

  input Modelica.SIunits.Time Ts = 0.1 "Sample time";
  input Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact "Discretization method";

  input Integer iu=1 "index of inout";
  input Integer iy=1 "index of output";
  output Boolean ok;

protected
  StateSpace ss=if systemOnFile then Modelica_LinearSystems2.StateSpace.Import.fromFile(fileName) else
      StateSpace(A=A, B=B, C=C, D=D);

 Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace dss=
                        Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace(
                                           ss,Ts,method);

algorithm
  assert(iu <= size(ss.B, 2) and iu > 0, "index for input is " + String(iu) +
    " which is not in [1, " + String(size(ss.B, 2)) + "].");
  assert(iy <= size(ss.C, 1) and iy > 0, "index for output is " + String(iy) +
    " which is not in [1, " + String(size(ss.C, 1)) + "].");

  Modelica_LinearSystems2.StateSpace.Plot.bodeSISO(
    ss,
    iu,
    iy);
  Modelica_LinearSystems2.WorkInProgress.DiscreteStateSpace.Plot.bodeSISO(
    dss,
    iu,
    iy);
  ok := true;

end plotBodeSISODiscrete;
