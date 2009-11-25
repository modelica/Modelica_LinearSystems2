within Modelica_LinearSystems2.Tests.Internal;
function getMatFiles "Returns a vector of files ending with .mat"
  extends Modelica.Icons.Function;
  input String directoryName = Modelica_LinearSystems2.Tests.TestDataDir annotation(Dialog);
  output String matFiles[:];
  output Integer nrMat;
protected
  Integer nEntries=Modelica.Utilities.Internal.FileSystem.getNumberOfFiles(directoryName);
  String files[nEntries];
  String files2[nEntries];
  Integer nrDir=nEntries;
  Integer i;
algorithm
  files := Modelica.Utilities.Internal.FileSystem.readDirectory(directoryName, nEntries);
  // count mat-files
  nrMat := 0;
  for i in 1:nrDir loop
    if Modelica.Utilities.Strings.find(files[i], ".mat", 1, false)<>0 then
      nrMat := nrMat + 1;
      files2[nrMat] := files[i];
    end if;
  end for;
  matFiles := files2[1:nrMat];

end getMatFiles;
