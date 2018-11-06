within Modelica_LinearSystems2.WorkInProgress.Tests.Analysis;
function fullAnalysis1
  "Example to compute the full analysis of a simple drive system"
  output Boolean ok;
algorithm
  Modelica_LinearSystems2.ModelAnalysis.FullAnalysis(
     "Modelica_LinearSystems2.WorkInProgress.Tests.Examples.SimpleDrive_SISO");
  ok := true;

end fullAnalysis1;
