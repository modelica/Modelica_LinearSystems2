within Modelica_LinearSystems2.Utilities.Plot;
function plotFFTs_fromDirectory "Plot amplitudes of FFT results (from result files in existing directory)"
  extends Modelica.Icons.Function;

  import Modelica.Utilities.Internal.FileSystem;
  import Modelica.Utilities.Streams.print;
  import Modelica.Utilities.Strings;
  import Modelica_LinearSystems2.Utilities.Plot;
  input String directory "Existing directory in which result data is present";
  input Boolean logX = false "= true, if logarithmic scale of x-axis" annotation(choices(checkBox=true));
  input Boolean fullPathTitle = false
    "= true, if directory path should be contained in the plot title, otherwise just file name"
    annotation(choices(checkBox=true));
protected
  String fftMarker = "FFT." "Substring at the beginning indicating a file containing fft";
  Integer nFiles = FileSystem.getNumberOfFiles(directory) "Number of files in the directory";
  String files[nFiles] = FileSystem.readDirectory(directory,nFiles) "List of files in the directory";
  String fft_files[nFiles] "List of files matching fftMarker condition";
  String fft_filesSorted[:];
  Integer nFFT=0 "Number of files matching fftMarker condition";
  String file;
  Real ix=0;
  Real iy=0;
  Real increment=10;
algorithm
  // Determine FFT files in directory
  for i in 1:nFiles loop
    // Determine whether file starts with fftMarker
    if Strings.length(files[i]) > Strings.length(fftMarker) then
      if Strings.isEqual(Strings.substring(files[i],1,Strings.length(fftMarker)), fftMarker) then
        // Store file in "fft_files"
        nFFT :=nFFT + 1;
        fft_files[nFFT] :=files[i];
      end if;
    end if;
  end for;

  if nFFT > 0 then
    // Sort the files
    fft_filesSorted := Modelica.Utilities.Strings.sort(fft_files[1:nFFT]);

    // Plot the files
    for i in 1:nFFT loop
      file := directory + "/" + fft_filesSorted[i];
      Modelica_LinearSystems2.Utilities.Plot.plotFFT_fromFile(
        file, logX, ix, iy, fullPathTitle);
      ix :=ix + increment;
      iy :=iy + increment;
    end for;
    print("Finished.");
  else
    print("No file for a FFT plot found in directory " + directory);
  end if;

  annotation (
    Documentation(revisions="<html>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>Date</th> <th align=\"left\">Description</th></tr>

<tr><td valign=\"top\"> Nov. 29, 2015 </td>
    <td valign=\"top\">
     Initial version implemented by
     Martin R. Kuhn and Martin Otter 
     (<a href=\"http://www.dlr.de/rmc/sr/en\">DLR Institute of System Dynamics and Control</a>)<br>  
     The research leading to these results has received funding from the European Union’s Seventh
     Framework Programme (FP7/2007-2016) for the Clean Sky Joint Technology Initiative under
     grant agreement no. CSJU-GAM-SGO-2008-001.</td></tr>
</table>
</html>", info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Utilities.Plot.plotFFTs_fromDirectory(
  directory, logX, fullPathTitle)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate plots of a Fast Fourier Transformation results of <b>all concerning files</b> saved in <code>directory</code>.
The files which contain the FFT results muss all start with a substring &quot;FFT.&quot;.
Only such a files will be proceeded by this function.
To generate the FFT result file, see e.g.
<a href=\"modelica://Modelica.Math.FastFourierTransform.Examples.RealFFT1\">Modelica.Math.FastFourierTransform.Examples.RealFFT1</a>.
</p>
</html>"));
end plotFFTs_fromDirectory;
