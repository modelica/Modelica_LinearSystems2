within Modelica_LinearSystems2.Internal;
record AnalyseOptions
  "Defines the characteristics of the eigenvalues to be print or to be plot"
  extends Modelica.Icons.Record;
  Boolean plotEigenValues=true "Plot eigenvalues"     annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotInvariantZeros=true "Plot invariant zeros"    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotStepResponse=true "Plot step respones: Only for SISO system"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotFrequencyResponse=true "Plot bode diagram: Only for SISO system"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printEigenValues=true "Write eigenvalues into the report"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printEigenValueProperties=true "Write eigenvalues with properties"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printInvariantZeros=true "Write invariant zreos into the report"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printControllability=true
    "Indicates controllability of every single pole"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printObservability=true
    "Indicates observability of every single pole"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  String headingEigenValues="Eigenvalues";
  String headingInvariantzeros="Invariant zeros";
  String headingStepResponse="Step response";
  String headingFrequencyResponse="Frequency response";

end AnalyseOptions;
