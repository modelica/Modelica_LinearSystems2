within Modelica_LinearSystems2.Internal.Streams;
record AnalyseOptions
  "Defines the characteristics of the eigenvalues to be print or to be plot"
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

end AnalyseOptions;
