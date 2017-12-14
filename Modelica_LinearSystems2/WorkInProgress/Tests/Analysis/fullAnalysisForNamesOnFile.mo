within Modelica_LinearSystems2.WorkInProgress.Tests.Analysis;
function fullAnalysisForNamesOnFile
  "Apply FullAnalysis on every class that is defined with its class name in a file"
  input String fileName=
      "modelica://Modelica_LinearSystems2.WorkInProgress.Tests.Analysis/blocksAndExamples.txt"
    "File with block and example names to be inspected";
  output Boolean ok;
protected
  String file[:]=Modelica.Utilities.Streams.readFile(Dymola_ResolveURI(fileName));
algorithm
  Modelica.Utilities.Files.removeFile("log.txt");
  for name in file loop
    Modelica.Utilities.Streams.print("fullAnalysis of " + name, "log.txt");
    Modelica_LinearSystems2.ModelAnalysis.FullAnalysis(name);
  end for;
  ok := true;

  annotation(__Dymola_interactive=true);
end fullAnalysisForNamesOnFile;
