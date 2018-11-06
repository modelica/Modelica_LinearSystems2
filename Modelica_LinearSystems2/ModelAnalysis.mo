within Modelica_LinearSystems2;
package ModelAnalysis
  "Package of functions to perform analysis on Modelica models (nonlinear models are linearized)"
  function Linearize
    "Linearize a model and return the linearized model as StateSpace object"
    extends Modelica_LinearSystems2.Internal.PartialAnalyzeFunction;
  public
    output Modelica_LinearSystems2.StateSpace ss = ssLin
      "Linearized system as StateSpace object";
  algorithm

    annotation (
      Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
            textString="L")}));
  end Linearize;
  extends Modelica.Icons.ExamplesPackage;
  function Poles "Linearize a model and plot the poles of the linearized model"
    extends Modelica_LinearSystems2.Internal.PartialAnalyzeFunction;
  algorithm
    Modelica_LinearSystems2.StateSpace.Plot.polesAndZeros(ssLin, zeros=false, print=true);
    annotation (
      Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
            textString=
                 "P")}));
  end Poles;

  function PolesAndZeros
    "Linearize a model and plot the poles and zeros of the linearized model"
    extends Modelica_LinearSystems2.Internal.PartialAnalyzeFunction;
  algorithm
    Modelica_LinearSystems2.StateSpace.Plot.polesAndZeros(ssLin, print=true);
    annotation (
      Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
            textString="PZ")}));
  end PolesAndZeros;

  function TransferFunctions
    "Linearize a model and plot the transfer functions from all inputs to all outputs of the linearized model"
    extends Modelica_LinearSystems2.Internal.PartialAnalyzeFunction;
    input Boolean dB_w = false
      "= true, to plot dB over w [rad/s] otherwise magnitude over f [Hz]" annotation(choices(checkBox=true));
  algorithm
    Modelica_LinearSystems2.StateSpace.Plot.bodeMIMO(ssLin,Hz=not dB_w, dB=dB_w);
    annotation (
      Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
            textString="TF")}));
  end TransferFunctions;

  function FullAnalysis
    "Linearize a model and perform all available linear analysis operations"
    extends Modelica_LinearSystems2.Internal.PartialAnalyzeFunction;
    input Boolean dB_w = false
      "= true, to plot Bode as dB over w [rad/s] otherwise magnitude over f [Hz]" annotation(choices(checkBox=true));
  algorithm
     Modelica_LinearSystems2.StateSpace.Analysis.analysis(
       ssLin,
       analyseOptions=Modelica_LinearSystems2.Internal.AnalyseOptions(dB_w=dB_w));
    annotation (
      Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
            textString="A")}));
  end FullAnalysis;

  function RootLocus = Modelica_LinearSystems2.Utilities.Plot.rootLocusOfModel
    "Compute and plot the root locus of one parameter of a model (= eigen values of the model that is linearized for every parameter value)"
    annotation (Icon(graphics={
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={255,127,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,-80}},
            lineColor={255,127,0},
          textString="RL")}));

end ModelAnalysis;
