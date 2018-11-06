within Modelica_LinearSystems2.Utilities.Import.Examples;
function linearizeDoublePendulum "Linearize double pendulum"
  extends Modelica.Icons.Function;

  output Real A[:,:] "A-matrix";
  output Real B[:,:] "B-matrix";
  output Real C[:,:] "C-matrix";
  output Real D[:,:] "D-matrix";
  output String inputNames[:] "Modelica names of inputs";
  output String outputNames[:] "Modelica names of outputs";
  output String stateNames[:] "Modelica names of states";
algorithm
  (A,B,C,D,inputNames,outputNames,stateNames) :=
    Modelica_LinearSystems2.Utilities.Import.linearize(
    "Modelica_LinearSystems2.Utilities.Plants.DoublePendulum",
    1.0);
  annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(A,B,C,D,inputNames,outputNames,stateNames) = Utilities.Import.Examples.linearizeDoublePendulum()
</pre></blockquote>

<h4>Description</h4>
<p>
Linearize <a href=\"Modelica_LinearSystems2.Utilities.Plants.DoublePendulum\">double pendulum model</a>.
The linearization is done at 1&nbsp;s of the simulation time of model.
</p>
</html>"));
end linearizeDoublePendulum;
