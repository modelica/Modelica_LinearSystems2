within Modelica_LinearSystems2.Internal;
partial function timeResponseMask_zp
  "Declares the common structure for the set of response functions for a zeros and poles transfer function"
  input ZerosAndPoles zp;
  input Real dt=0 "Sample time [s]";
  input Real tSpan=0 "Simulation time span [s]";
  input Integer toWindow[:,:]=fill(
      -1,
      1,
      1) "-1/0/>0 plot in new window/last window/window ID";
  /* the Dialog annotation seems not to work properly. why? */
  input Boolean clearWindow[:,:]=fill(
      false,
      1,
      1) "= true, if previous window content is removed"             annotation(Dialog(enable=max(toWindow)>=0));
  input String heading[:,:]=fill(
      "Time response",
      1,
      1) "Heading of the response diagram";

  input String columnLabels[2]=
      Modelica_LinearSystems2.Internal.defaultColumnLabels(1);

  replaceable output Real y[:,1,1]
    "Output response: (number of samples) x (number of outputs) x (number of inputs)";
  output Real t[:] "Time vector: (number of samples)";
  replaceable output Real x_continuous[:, Modelica_LinearSystems2.ZerosAndPoles.Analysis.denominatorDegree(zp), 1]
    "State trajectories: (number of samples) x (number of states) x (number of inputs)";

end timeResponseMask_zp;
