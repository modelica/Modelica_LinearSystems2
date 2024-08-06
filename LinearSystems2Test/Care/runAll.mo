within LinearSystems2Test.Care;
function runAll
  "Executes all functions with name care*.mo in this directory/package"
  extends Modelica.Icons.Function;

  input String logFile = "";
  output Boolean ok;

algorithm
  ok := LinearSystems2Test.Care.care1(logFile);
  ok := LinearSystems2Test.Care.care2(logFile);
  ok := LinearSystems2Test.Care.care3(logFile);
  ok := LinearSystems2Test.Care.care4(logFile);
  ok := LinearSystems2Test.Care.care5(logFile);
  ok := LinearSystems2Test.Care.care6(logFile);
  ok := LinearSystems2Test.Care.care7(logFile);
  ok := LinearSystems2Test.Care.care8(logFile);
  ok := LinearSystems2Test.Care.care9(logFile);
  ok := LinearSystems2Test.Care.care10(logFile);
  ok := LinearSystems2Test.Care.care11(logFile);
  ok := LinearSystems2Test.Care.care12(logFile);
  ok := LinearSystems2Test.Care.care13(logFile);
  ok := LinearSystems2Test.Care.care14(logFile);
  ok := LinearSystems2Test.Care.care15(logFile);
end runAll;
