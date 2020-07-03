within Modelica_LinearSystems2.Utilities.Plot;
function plot_FFT_fromFile
  "Plot amplitudes of FFT results (from result file)"
  import Modelica_LinearSystems2.Utilities.Plot;
  input String fileName "File where FFT data is stored"
    annotation (Dialog(
      loadSelector(filter="MATLAB MAT-files (*.mat)",
          caption="Open file of FFT result data")));
  input Boolean logX = false "= true, if logarithmic scale of x-axis" annotation(choices(checkBox=true));
protected
   Integer dims[2] = Modelica.Utilities.Streams.readMatrixSize(
     fileName,"FFT");
   Real fA[:,:] = Modelica.Utilities.Streams.readRealMatrix(
     fileName, "FFT", dims[1], dims[2]);
algorithm
  Plot.diagram(Plot.Records.Diagram(curve={Plot.Records.Curve(x=if logX then fA[4:size(fA,1),1] else fA[:,1], y=if logX then fA[4:size(fA,1),2] else fA[:,2])},
                                    heading="Result of FFT calculation (" + fileName+")",
                                    xLabel="Frequency in [Hz]", yLabel="Amplitude", logX=logX));
  annotation(__Dymola_interactive=true, Documentation(revisions="<html>
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
</html>"));
end plot_FFT_fromFile;
