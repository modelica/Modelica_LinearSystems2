within Modelica_LinearSystems2.Internal;
partial function timeResponseMask_discrete
  "Declares the common structure for the set of response functions"
  extends Modelica.Icons.Function;

  input DiscreteStateSpace dss "Discrete state space system";
  input Real tSpan=0 "Simulation time span [s]";

  replaceable output Real y[:,size(dss.C, 1),size(dss.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_discrete[:,size(dss.A, 1),size(dss.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask_discrete;
