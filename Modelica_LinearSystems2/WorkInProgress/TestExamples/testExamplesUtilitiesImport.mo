within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesUtilitiesImport "Test all examples from package Utilities.Import"

algorithm
  Modelica_LinearSystems2.Utilities.Import.Examples.linearizeDoublePendulum();

  Modelica.Utilities.Streams.print("testExamplesUtilitiesImport done!");
end testExamplesUtilitiesImport;
