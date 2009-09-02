within Modelica_LinearSystems2.Internal;
partial function timeResponseMask
  "Declares the common structure for the set of response functions"
  input StateSpace sc;
  input Real dt=0 "Sample time [s]";
  input Real tSpan=0 "Simulation time span [s]";
  input Integer toWindow[:,:]=fill(
      -1,
      size(sc.C, 1),
      size(sc.B, 2)) "-1/0/>0 plot in new window/last window/window ID";
  /* the Dialog annotation seems not to work properly. why? */
  input Boolean clearWindow[:,:]=fill(
      false,
      size(sc.C, 1),
      size(sc.B, 2)) "= true, if previous window content is removed" annotation(Dialog(enable=max(toWindow)>=0));
  input String heading[:,:]=fill(
      "Time response",
      size(sc.C, 1),
      size(sc.B, 2)) "Heading of the response diagram";

  input String columnLabels[size(sc.C, 1) + 1]=
      Modelica_LinearSystems2.Internal.defaultColumnLabels(size(sc.C, 1));

  replaceable output Real y[:,size(sc.C, 1),size(sc.B, 2)]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_continuous[:,size(sc.A, 1),size(sc.B, 2)]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask;
