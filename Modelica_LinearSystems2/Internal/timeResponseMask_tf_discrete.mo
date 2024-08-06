within Modelica_LinearSystems2.Internal;
partial function timeResponseMask_tf_discrete
  "Declares the common structure for the set of discrete response functions"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.DiscreteTransferFunction dtf;
  input Real tSpan=0 "Simulation time span [s]";

  replaceable output Real y[:,1,1]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_discrete[:,Modelica_LinearSystems2.DiscreteTransferFunction.Analysis.denominatorDegree(dtf), 1]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask_tf_discrete;
