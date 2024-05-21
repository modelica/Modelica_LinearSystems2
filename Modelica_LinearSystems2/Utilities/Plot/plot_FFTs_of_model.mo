within Modelica_LinearSystems2.Utilities.Plot;
function plot_FFTs_of_model
  "Plot amplitudes of FFT results (from result files of translated model)"
  extends Modelica.Icons.Function;

  input String modelName
    "Model that was used to generate FFT data (FFT result files are stored in directory <modelName>)" annotation(Dialog(__Dymola_translatedModel=true));
  input Boolean logX = false "= trrue, if logarithmic scale of x-axis" annotation(choices(checkBox=true));
protected
  String directory = Internal.getLastName(modelName);

algorithm
  // Check if directory exists
  if not Modelica.Utilities.Files.exist(directory) then
    Modelica.Utilities.Streams.print(
      "... Do not find FFT result files, because no directory : " +
      Modelica.Utilities.Files.fullPathName(directory));
    return;
  end if;

  // Read FFT.* files and plot them
  plot_FFTs_from_directory(directory=directory, logX=logX);

  annotation(__Dymola_interactive=true, Documentation(revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
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
end plot_FFTs_of_model;
