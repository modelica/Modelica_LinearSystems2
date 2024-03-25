within Modelica_LinearSystems2.Internal;
record AnalyseOptions
  "Defines the characteristics of the eigenvalues to be print or to be plot"
  extends Modelica.Icons.Record;
  Boolean plotEigenValues=true "Plot eigenvalues" annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotInvariantZeros=true "Plot invariant zeros" annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotStepResponse=true "Plot step response"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean plotFrequencyResponse=true "Plot bode diagram"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));
  Boolean printSystem=true "Write system into the report (if not too large)"
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
  Boolean dB_w = false
    "= true, to plot Bode with dB over w [rad/s] otherwise magnitude over f [Hz]"
    annotation(Dialog(group="Analyse options"),choices(checkBox=true));

end AnalyseOptions;
