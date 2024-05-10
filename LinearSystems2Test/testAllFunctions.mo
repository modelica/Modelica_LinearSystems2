within LinearSystems2Test;
function testAllFunctions "Runs all test cases for functions"
  extends Modelica.Icons.Function;
  import Modelica.Utilities.Streams.print;
  input String logFile = "LinearSystems2TestLog.txt"
    "Filename where the log of all functions is stored";
  output Boolean ok;
protected
  Boolean result;
  String file;
algorithm
  ok:=false;
  file :=Modelica.Utilities.Files.fullPathName(logFile);
  print("... testAllFunctions(..) is logged in " + file);

  if not Modelica.Utilities.Strings.isEmpty(file) then
    Modelica.Utilities.Files.removeFile(file);
  end if;

  print("--- Test functions of LinearSystems2 library");
  print("--- Test functions of LinearSystems2 library", logFile);
  print("--- Test functions of Modelica_LinearSystems2.Examples.StateSpace");
  Modelica_LinearSystems2.Examples.StateSpace.analysis();
  Modelica_LinearSystems2.Examples.StateSpace.analysisDcGain();
  Modelica_LinearSystems2.Examples.StateSpace.analysisControllability();
  Modelica_LinearSystems2.Examples.StateSpace.analysisObservability();
  Modelica_LinearSystems2.Examples.StateSpace.analysisControllablePoles();
  Modelica_LinearSystems2.Examples.StateSpace.analysisTimeResponse();
  Modelica_LinearSystems2.Examples.StateSpace.analysisInitialResponse();
  Modelica_LinearSystems2.Examples.StateSpace.analysisImpulseResponse();
  Modelica_LinearSystems2.Examples.StateSpace.analysisInvariantZeros();
  Modelica_LinearSystems2.Examples.StateSpace.analysisPolesAndZerosSISO();
  Modelica_LinearSystems2.Examples.StateSpace.analysisStepResponse();
  Modelica_LinearSystems2.Examples.StateSpace.analysisStairCase();
  Modelica_LinearSystems2.Examples.StateSpace.conversionFromTransferFunction();
  Modelica_LinearSystems2.Examples.StateSpace.conversionFromZerosAndPoles();
  Modelica_LinearSystems2.Examples.StateSpace.conversionToTransferFunctionMIMO();
  Modelica_LinearSystems2.Examples.StateSpace.conversionToTransferFunctionSISO();
  Modelica_LinearSystems2.Examples.StateSpace.conversionToZerosAndPolesMIMO();
  Modelica_LinearSystems2.Examples.StateSpace.conversionToZerosAndPolesSISO();
  Modelica_LinearSystems2.Examples.StateSpace.designAssignPolesSISO();
  Modelica_LinearSystems2.Examples.StateSpace.designAssignPolesMIMO();
  Modelica_LinearSystems2.Examples.StateSpace.designKalmanFilter();
  Modelica_LinearSystems2.Examples.StateSpace.designLQG();
  Modelica_LinearSystems2.Examples.StateSpace.designLQR();
  Modelica_LinearSystems2.Examples.StateSpace.importFromModel();
  Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros();
  Modelica_LinearSystems2.Examples.StateSpace.plotPolesAndZeros2();
  Modelica_LinearSystems2.Examples.StateSpace.plotBodeMIMO();
  Modelica_LinearSystems2.Examples.StateSpace.plotBodeSISO();
  Modelica_LinearSystems2.Examples.StateSpace.plotImpulse();
  Modelica_LinearSystems2.Examples.StateSpace.plotInital();
  Modelica_LinearSystems2.Examples.StateSpace.plotRamp();
  Modelica_LinearSystems2.Examples.StateSpace.plotStep();
  Modelica_LinearSystems2.Examples.StateSpace.plotTimeResponse();
  Modelica_LinearSystems2.Examples.StateSpace.plotZeros();
  Modelica_LinearSystems2.Examples.StateSpace.transformation();
  Modelica_LinearSystems2.Examples.StateSpace.transformationExtract();
  Modelica_LinearSystems2.Examples.StateSpace.transformationToIrreducibleForm();

  print("--- Test functions of Modelica_LinearSystems2.Examples.TransferFunction");
  Modelica_LinearSystems2.Examples.TransferFunction.plotPolesAndZeros();
  Modelica_LinearSystems2.Examples.TransferFunction.importFromFile();
  Modelica_LinearSystems2.Examples.TransferFunction.operationsOnTransferFunctions();
  Modelica_LinearSystems2.Examples.TransferFunction.plotBode1();
  Modelica_LinearSystems2.Examples.TransferFunction.plotBode2();
  Modelica_LinearSystems2.Examples.TransferFunction.plotBode3();
  Modelica_LinearSystems2.Examples.TransferFunction.plotBodeFilter2();
  Modelica_LinearSystems2.Examples.TransferFunction.plotImpulse();
  Modelica_LinearSystems2.Examples.TransferFunction.plotInital();
  Modelica_LinearSystems2.Examples.TransferFunction.plotStep();

  print("--- Test functions of Modelica_LinearSystems2.Examples.ZerosAndPoles");
  Modelica_LinearSystems2.Examples.ZerosAndPoles.analysisZerosAndPoles();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.analysisDcGain();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.conversionToStateSpace();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotPolesAndZeros();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBode1();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBode2();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBode3();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter1();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter2();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotBodeFilter3();
  Modelica_LinearSystems2.Examples.ZerosAndPoles.plotStep();

  print("--- Test functions of Modelica_LinearSystems2.Examples.DiscreteStateSpace");
  Modelica_LinearSystems2.Examples.DiscreteStateSpace.analysisEigenvalues();
  Modelica_LinearSystems2.Examples.DiscreteStateSpace.analysisTimeResponse();
  Modelica_LinearSystems2.Examples.DiscreteStateSpace.importFromModel();
  Modelica_LinearSystems2.Examples.DiscreteStateSpace.plotBodeSISO();

  print("--- Test functions of Modelica_LinearSystems2.Examples.DiscreteTransferFunction");
  Modelica_LinearSystems2.Examples.DiscreteTransferFunction.plotBode();

  print("--- Test functions of Modelica_LinearSystems2.Examples.DiscreteZerosAndPoles");
  Modelica_LinearSystems2.Examples.DiscreteZerosAndPoles.plotBode();

//   print("--- Test functions of LinearSystems2Test library");
//   print("--- Test functions of LinearSystems2Test library", logFile);
//   result := ModelicaTest.Math.ScalarFunctions(logFile);
  print("--- Test functions run done");

  ok := true;
end testAllFunctions;
