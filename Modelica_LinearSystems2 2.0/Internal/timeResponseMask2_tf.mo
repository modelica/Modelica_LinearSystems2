within Modelica_LinearSystems2.Internal;
partial function timeResponseMask2_tf
  "Declares the common structure for the set of response functions"
  input Modelica_LinearSystems2.TransferFunction tf;
  input Real dt=0 "Sample time [s]";
  input Real tSpan=0 "Simulation time span [s]";

  replaceable output Real y[:,1,1]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_continuous[:,Modelica_LinearSystems2.TransferFunction.Analysis.denominatorDegree(tf), 1]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask2_tf;
