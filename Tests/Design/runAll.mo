within Modelica_LinearSystems2.Tests.Design;
function runAll
  "Executes controller design for the according data in directory Test\\data"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Tests;

  input String dataSpecifier = "data_";
  input String fileSpecifier = "Design";
  input String directoryName= classDirectory()+ "../" + "Data";
  input String outputFile = fileSpecifier+".txt";
  input Types.AssignPolesMethod method=Tests.Types.AssignPolesMethod.KNV
    "method for pole assignment";
  output Integer numberOfFunctions;
  output Boolean ok;

protected
  String functions[:];
  String functions2[:];
  String designMos = fileSpecifier + ".mos";
  String h="Design.";
  Integer index;
  Integer n;
  Integer i;

algorithm
  if Modelica.Utilities.Files.exist(designMos) then
    Modelica.Utilities.Files.removeFile(designMos);
  end if;
  if Modelica.Utilities.Files.exist(outputFile) then
    Modelica.Utilities.Files.removeFile(outputFile);
  end if;

  Modelica.Utilities.Streams.print("import Modelica_LinearSystems2.Tests",designMos);
  Modelica.Utilities.Streams.print("import Modelica_LinearSystems2.Tests.Types",designMos);

   (functions, n) := Tests.Internal.getMatFiles(dataSpecifier,directoryName);
   functions2 := fill("",n);
   for i in 1:n loop
    index := Modelica.Utilities.Strings.find(functions[i],".");
    functions2[i] := Modelica.Utilities.Strings.substring(functions[i],1,index-1);
    functions[i]:="/"+functions[i];
    if method==Tests.Types.AssignPolesMethod.KNV then
      Modelica.Utilities.Streams.print("Modelica.Utilities.Streams.print(\"\\n\\n### "+functions2[i]+" KNV-method ###\\n\",\"" + outputFile+ "\")",designMos);
      Modelica.Utilities.Streams.print("Tests.Design.testPoleAssignment2(\""+directoryName+functions[i]+"\", Types.AssignPolesMethod.KNV, true ,\""+outputFile+"\", false);",designMos);
    else
      Modelica.Utilities.Streams.print("Modelica.Utilities.Streams.print(\"\\n\\n### "+functions2[i]+" Schur-method ###\\n\",\"" + outputFile+ "\")",designMos);
      Modelica.Utilities.Streams.print("Tests.Design.testPoleAssignment2(\""+directoryName+functions[i]+"\", Types.AssignPolesMethod.Schur, true,\""+outputFile+"\", false);",designMos);
    end if;
  end for;
   ok := RunScript(designMos);
//ok := true;

numberOfFunctions := n;
end runAll;
