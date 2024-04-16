within Modelica_LinearSystems2.Utilities.Plot;
function plot_FFTs_from_directory
  "Plot amplitudes of FFT results (from result files in existing directory)"
  import Modelica.Utilities.Internal.FileSystem;
  import Modelica.Utilities.Strings;
  import Modelica_LinearSystems2.Utilities.Plot;
  input String directory
    "Existing directory in which result data is present";
  input Boolean logX = false "= trrue, if logarithmic scale of x-axis" annotation(choices(checkBox=true));
protected
  Integer nFiles = FileSystem.getNumberOfFiles(directory);
  String files[nFiles] = FileSystem.readDirectory(directory,nFiles);
  String fft_files[nFiles];
  String fft_filesSorted[:];
  Integer nFFT=0;
  Integer dims[2];
  Real fA[:,:];
  String file;
  Real ix=0;
  Real iy=0;
  Real increment=10;
algorithm
  // Determine FFT files in directory
  for i in 1:nFiles loop
     // Determine whether file starts with "FFT."
     if Strings.length(files[i]) > 4 then
        if Strings.substring(files[i],1,4) == "FFT." then
           // Store file in "fft_files"
           nFFT :=nFFT + 1;
           fft_files[nFFT] :=files[i];
        end if;
     end if;
  end for;

  // Sort the files
  fft_filesSorted := Modelica.Utilities.Strings.sort(fft_files[1:nFFT]);

  // Plot the files
  for i in 1:nFFT loop
     file := directory + "/" + fft_filesSorted[i];
     dims := Modelica.Utilities.Streams.readMatrixSize(file, "FFT");
     fA   := Modelica.Utilities.Streams.readRealMatrix(file, "FFT", dims[1], dims[2]);

     ix :=ix + increment;
     iy :=iy + increment;
     Plot.diagram(Plot.Records.Diagram(curve={
          Plot.Records.Curve(x=if logX then fA[4:size(fA,1),1] else fA[:,1],
                             y=if logX then fA[4:size(fA,1),2] else fA[:,2])},
          heading="Result of FFT calculation (" + file +")",
          xLabel="Frequency in [Hz]", yLabel="Amplitude", logX=logX),
          Plot.Records.Device(xTopLeft=ix,
                              yTopLeft=iy));
  end for;
  annotation(__Dymola_interactive=true, Documentation(revisions="<html>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td valign=\"top\"> Nov. 29, 2015 </td>
    <td valign=\"top\">
     Initial version implemented by
     Martin R. Kuhn and Martin Otter 
     (<a href=\"https://www.dlr.de/sr/en/\">DLR Institute of System Dynamics and Control</a>)<br>  
     The research leading to these results has received funding from the European Union’s Seventh
     Framework Programme (FP7/2007-2016) for the Clean Sky Joint Technology Initiative under
     grant agreement no. CSJU-GAM-SGO-2008-001.</td></tr>
</table>
</html>"));
end plot_FFTs_from_directory;
