within Modelica_LinearSystems2.WorkInProgress.TestExamples;
function testExamplesMathMatrices "Test all examples from package Math.Matrices"

algorithm
  Modelica_LinearSystems2.Math.Matrices.Examples.exampleHessenberg();
  Modelica_LinearSystems2.Math.Matrices.Examples.exampleQR();
  Modelica_LinearSystems2.Math.Matrices.Examples.exampleSVD();
  Modelica_LinearSystems2.Math.Matrices.Examples.care();

  Modelica.Utilities.Streams.print("testExamplesMathMatrices done!");
end testExamplesMathMatrices;
