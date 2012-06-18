within Modelica_LinearSystems2;
package Controller "Continuous and discrete input/output blocks. Easy to switch from continuous to discrete representation."
  extends Modelica.Icons.Package;

constant String DataDir=Modelica_LinearSystems2.DataDir
  "Absolute path to directory containing utilitiy files for this package, such as images";


  extends Modelica.Icons.Package;


  annotation (Documentation(info="<html>
<p>
This library provides input/output blocks where every
block is available in a <b>continuous</b> and a <b>discrete</b> (sampled) 
representation. A block is defined via its <b>continuous
parameterization</b>. By specifying a discretization method and 
a sample time, the discrete representation is automatically
derived from the continuous form. The defaults of the most
important options for <b>all blocks</b> are set in the global SampleClock 
component (via inner/outer). 
As a result, it is, e.g., easy to switch quickly 
between a continuous and a discrete representation of all 
blocks of a controller. 
</p>
 
<p>
Examples to demonstrate the technique are given in sublibrary
<a href=\"modelica://Modelica_LinearSystems2.Controller.Examples\">Examples</a>.
Especially, the continuous or discrete control of a simple flexible
drive with a P-PI cascade controller is demonstrated in example
<a href=\"modelica://Modelica_LinearSystems2.Controller.Examples.SimpleControlledDrive\">SimpleControlledDrive</a>.
</p>
 
<p align=\"center\">
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/SimpleControlledDrive_Plot2.png\">
</p>
 
</html>"));
end Controller;
