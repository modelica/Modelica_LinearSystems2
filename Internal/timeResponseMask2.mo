within Modelica_LinearSystems2.Internal;
partial function timeResponseMask2
  "Declares the common structure for the set of response functions"
  input StateSpace sc;
  input Real dt=0 "Sample time [s]";
  input Real tSpan=0 "Simulation time span [s]";

  replaceable output Real y[:,size(sc.C, 1),size(sc.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_continuous[:,size(sc.A, 1),size(sc.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask2;
