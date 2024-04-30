within Modelica_LinearSystems2.Utilities.Streams;
record AnalyseOptions
  "Defines the characteristics of the eigenvalues for print or plot"
  extends Modelica.Icons.Record;
  Boolean plotEigenValues = true;
  Boolean plotInvariantZeros = true;
  String headingEigenValues = "Eigen values and invariant zeros";
  Boolean printEigenValues = true;
  Boolean printInvariantZeros = true;
  Boolean printControllability = true;
  Boolean printObservability = true;
  Boolean plotStepResponse = true "only for SISO system";
  String headingStepResponse = "Step response";
  Boolean plotFrequencyResponse = true "only for SISO system";
  String headingFrequencyResponse = "Frequency response";

  annotation (Documentation(info="<html>
<p>
Record collecting parameters for the characteristics of eigenvalues to be printed or plotted.
</p>
</html>"));
end AnalyseOptions;
