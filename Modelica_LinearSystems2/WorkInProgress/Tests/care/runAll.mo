within Modelica_LinearSystems2.WorkInProgress.Tests.care;
function runAll
  "Executes all functions with name care*.mo in this directory/package"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.WorkInProgress.Tests;

  input String specifier = "care";
  input String outputFile = specifier +".txt";
  output Integer numberOfFunctions;
  output Boolean ok;

protected
  String careMos = specifier + ".mos";
  String functions[:];
  String functions2[:];
  String h="care.";
  Integer index;
  Integer n;
  Integer i;

algorithm
  if Modelica.Utilities.Files.exist(careMos) then
    Modelica.Utilities.Files.removeFile(careMos);
  end if;
  if Modelica.Utilities.Files.exist(outputFile) then
    Modelica.Utilities.Files.removeFile(outputFile);
  end if;
  Modelica.Utilities.Streams.print("import Modelica_LinearSystems2.WorkInProgress.Tests.care",careMos);

  (functions, n) := Tests.Internal.getMoFiles(specifier);
  functions2 := fill("",n);
  for i in 1:n loop
    index := Modelica.Utilities.Strings.find(functions[i],".");
    functions2[i] := Modelica.Utilities.Strings.substring(functions[i],1,index-1);
    Modelica.Utilities.Streams.print("Modelica.Utilities.Streams.print(\"\\n\\n### "+functions2[i]+" ###\\n\",\"" + outputFile+ "\")",careMos);
    functions2[i] := h+functions2[i] +"(\""+outputFile+"\" )";
    Modelica.Utilities.Streams.print(functions2[i],careMos);
  end for;
    ok := DymolaCommands.SimulatorAPI.RunScript(careMos);

  numberOfFunctions := n;

end runAll;
