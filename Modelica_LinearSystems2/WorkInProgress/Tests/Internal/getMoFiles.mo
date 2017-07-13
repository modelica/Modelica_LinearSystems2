within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function getMoFiles "Returns a vector of files with *.mo"
  extends Modelica.Icons.Function;
  input String specifier="care" "Subdirectory and filename specifier";
  input String directoryName=
    Modelica.Utilities.Files.loadResource("modelica://Modelica_LinearSystems2/WorkInProgress/Tests/" + specifier)
    annotation(Dialog);
  output String moFiles[:];
  output Integer nrMo;
protected
  Integer nEntries=Modelica.Utilities.Internal.FileSystem.getNumberOfFiles(directoryName);
  String files[nEntries];
  String files2[nEntries];
  Integer nrDir=nEntries;
  Integer i;
algorithm
  files := Modelica.Utilities.Internal.FileSystem.readDirectory(directoryName, nEntries);
  // count mo-files
  nrMo := 0;
  for i in 1:nrDir loop
    if Modelica.Utilities.Strings.find(files[i], specifier, 1, false)<>0 and Modelica.Utilities.Strings.find(files[i], ".mo", 1, false)<>0 then
      nrMo := nrMo + 1;
      files2[nrMo] := files[i];
    end if;
  end for;
  moFiles := files2[1:nrMo];

end getMoFiles;
