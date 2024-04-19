within Modelica_LinearSystems2;
package Controllers "Package of continuous and discrete input/output blocks"
  extends Modelica.Icons.Package;

  constant String DataDir=Modelica_LinearSystems2.DataDir
    "Absolute path to directory containing utility files for this package, such as images";


  annotation (
    Documentation(
      info="<html>
<p>
This package provides input/output blocks where every
block is available in a <strong>continuous</strong> and a <strong>discrete</strong> (sampled)
representation. A block is defined via its <strong>continuous
parameterization</strong>. By specifying a discretization method and
a sample time, the discrete representation is automatically
derived from the continuous form. The defaults of the most
important options for <strong>all blocks</strong> are set in the global SampleClock
component (via inner/outer).
As a result, it is, e.g., easy to switch quickly
between a continuous and a discrete representation of all
blocks of a controller.
</p>

<p>
Examples to demonstrate the technique are given in the package
<a href=\"modelica://Modelica_LinearSystems2.Controllers.Examples\">Examples</a>.
Especially, the continuous or discrete control of a simple flexible
drive with a P-PI cascade controller is demonstrated in example
<a href=\"modelica://Modelica_LinearSystems2.Controllers.Examples.SimpleControlledDrive\">SimpleControlledDrive</a>.
</p>

<p>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/Examples/SimpleControlledDrive_Layout.png\">
</p>
</html>"));
end Controllers;
