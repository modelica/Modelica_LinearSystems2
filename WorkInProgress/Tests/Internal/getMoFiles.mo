within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function getMoFiles "Returns a vector of files with *.mo"
  extends Modelica.Icons.Function;
  input String specifier="data";
  input String directoryName = classDirectory()+ "../" + specifier annotation(Dialog);
  output String moFiles[:];
  output Integer nrMat;
protected
  Integer nEntries=Modelica.Utilities.Internal.FileSystem.getNumberOfFiles(directoryName);
  String files[nEntries];
  String files2[nEntries];
  Integer nrDir=nEntries;
  Integer i;
algorithm
  files := Modelica.Utilities.Internal.FileSystem.readDirectory(directoryName, nEntries);
  // count mo-files
  nrMat := 0;
  for i in 1:nrDir loop
    if Modelica.Utilities.Strings.find(files[i], specifier, 1, false)<>0 and Modelica.Utilities.Strings.find(files[i], ".mo", 1, false)<>0 then
      nrMat := nrMat + 1;
      files2[nrMat] := files[i];
    end if;
  end for;
  moFiles := files2[1:nrMat];

end getMoFiles;
