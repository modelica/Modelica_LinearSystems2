within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesUtilitiesPlot "Test all examples from package Utilities.Plot"

algorithm
  Modelica_LinearSystems2.Utilities.Plot.Examples.plotSine();
  Modelica_LinearSystems2.Utilities.Plot.Examples.plotTwoSine();
  Modelica_LinearSystems2.Utilities.Plot.Examples.plotTwoSineDifferentStyles();
  Modelica_LinearSystems2.Utilities.Plot.Examples.showSinesInVectorDiagrams();
  Modelica_LinearSystems2.Utilities.Plot.Examples.showMatrixDiagrams();
  Modelica_LinearSystems2.Utilities.Plot.Examples.plotParameterizedCurve1();
  Modelica_LinearSystems2.Utilities.Plot.Examples.plotParameterizedCurve2();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfDrive();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfDriveSolidLine();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfControlledSISO1Log();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfControlledSISO2();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfControlledSISO2Log();
  Modelica_LinearSystems2.Utilities.Plot.Examples.rootLocusOfPIDDrive();
  Modelica_LinearSystems2.Utilities.Plot.Examples.showLinePatterns();
  Modelica_LinearSystems2.Utilities.Plot.Examples.showLegendStyles();
  Modelica_LinearSystems2.Utilities.Plot.Examples.showPointSymbols();

  Modelica.Utilities.Streams.print("testExamplesUtilitiesPlot done!");
end testExamplesUtilitiesPlot;
