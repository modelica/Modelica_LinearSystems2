within Modelica_LinearSystems2.Utilities.Plot;
function plotFFT_fromFile "Plot amplitudes of FFT results (from result file)"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Utilities.Plot;
  input String fileName "File where FFT data is stored"
    annotation (Dialog(
      loadSelector(
        filter="MATLAB MAT-files (*.mat)",
        caption="Open file of FFT result data")));
  input Boolean logX = false "= true, if logarithmic scale of x-axis" annotation(choices(checkBox=true));
  input Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm xTopLeft=0
    "Horizontal position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm yTopLeft=0
    "Vertical position of top left figure corner if applicable (e.g. window)"
    annotation(Dialog);
  input Boolean fullPathTitle = false
    "= true, if directory path should be contained in the plot title, otherwise just file name"
    annotation(choices(checkBox=true));
protected
  Integer dims[2] = Modelica.Utilities.Streams.readMatrixSize(fileName,"FFT")
    "Dimension of FFT matrix in file";
  Real fA[:,:] = Modelica.Utilities.Streams.readRealMatrix(fileName, "FFT", dims[1], dims[2])
    "Frequency and amplitude";
  String splitPath[3] "fileName splitted into directory, file name and extension";
  String fileNamePrint "Name of file to be contained in plot title";
algorithm
  // How file name to be printed in the plot title
  (splitPath[1],splitPath[2],splitPath[3]) := Modelica.Utilities.Files.splitPathName(fileName);
  if fullPathTitle then
    fileNamePrint := Modelica.Utilities.Files.fullPathName(fileName);
  else
    fileNamePrint := splitPath[2] + splitPath[3];
  end if;

  Modelica.Utilities.Streams.print("... generating FFT plot for " + fileName);

  Plot.diagram(
    Plot.Records.Diagram(
      curve={
        Plot.Records.Curve(
          x=if logX then fA[4:size(fA,1),1] else fA[:,1],
          y=if logX then fA[4:size(fA,1),2] else fA[:,2])},
      heading="Result of FFT calculation (" + fileNamePrint+")",
      xLabel="Frequency [Hz]",
      yLabel="Amplitude",
      logX=logX),
    Plot.Records.Device(
      xTopLeft=xTopLeft,
      yTopLeft=yTopLeft));

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
Utilities.Plot.plotFFT_fromFile(
  fileName,logX,xTopLeft,yTopLeft,fullPathTitle)
</pre></blockquote>

<h4>Description</h4>
<p>
Generate plot of a Fast Fourier Transformation results saved on a file <code>fileName</code>.
To generate the FFT result file, see e.g.
<a href=\"modelica://Modelica.Math.FastFourierTransform.Examples.RealFFT1\">Modelica.Math.FastFourierTransform.Examples.RealFFT1</a>.
</p>
</html>"));
end plotFFT_fromFile;
