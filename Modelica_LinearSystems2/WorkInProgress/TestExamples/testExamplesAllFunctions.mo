within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesAllFunctions "Call all functions across the library defined as examples"
  extends Modelica.Icons.Function;
algorithm
  Modelica_LinearSystems2.WorkInProgress.TestExamples.testExamplesMain();
  Modelica_LinearSystems2.WorkInProgress.TestExamples.testExamplesMathMatrices();
  Modelica_LinearSystems2.WorkInProgress.TestExamples.testExamplesUtilitiesPlot();
  Modelica_LinearSystems2.WorkInProgress.TestExamples.testExamplesUtilitiesImport();
  //  Modelica_LinearSystems2.WorkInProgress.TestExamples.testExamplesWorkInProgress();
end testExamplesAllFunctions;
