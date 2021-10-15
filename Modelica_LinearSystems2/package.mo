within ;
package Modelica_LinearSystems2 "Modelica_LinearSystems2 (version 2.4.1-beta.2) - Analysis, Synthesis and Modeling of Continuous and Discrete Linear Systems"

  extends Modelica.Icons.Package;

  constant String DataDir=classDirectory() +  "Resources/Data/"
  "Absolute path to directory containing utilitiy files for this package";


annotation (
  preferredView="info",
  uses(
    Modelica(version="4.0.0"),
    DymolaCommands(version="1.11")),
  version="2.4.1-beta.2",
  versionDate="2021-10-29",
  dateModified = "2021-10-15 14:00:00Z",
  revisionId="$F​ormat:%h %ci$",
  conversion(
    from(version={"2.0", "2.1", "2.2", "2.3", "2.3.1", "2.3.2", "2.3.2", "2.3.3", "2.3.4"},
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.4.mos", to="2.3.5"),
    from(version="2.3.5",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.5.mos"),
    from(version="2.4.0",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.4.0.mos")),
  Documentation(info="<html>
<p>
Library <b>Modelica_LinearSystems2</b> is a Modelica package
providing different representations of linear, time invariant differential and
difference equation systems. For example, record
<a href=\"modelica://Modelica_LinearSystems2.StateSpace\">StateSpace</a>
defines a linear time invariant differential
equation system in state space form:
</p>
<blockquote><pre>
der(<b>x</b>) = <b>A</b> * <b>x</b> + <b>B</b> * <b>u</b>
    <b>y</b>  = <b>C</b> * <b>x</b> + <b>D</b> * <b>u</b>
</pre></blockquote>
<p>
Operators are overloaded to work conveniently with these system descriptions in an
interactive environment, e.g. to multiply transfer functions or to operate on complex numbers.
About 180 functions are provided to operate
on these data structures, e.g. to compute eigenvalues, zeros, step responses,
to design pole-placement and LQG controllers, to plot step responses, frequency responses,
eigenvalues, to convert between different description forms, or to
generate a linear system description by linearization of a Modelica model.
</p>

<p>
Furthermore, in subpackage
<a href=\"modelica://Modelica_LinearSystems2.Controller\">Controller</a>
about 20 input/output blocks of linear systems are provided that are
based on the different representation forms, e.g. PID, StateSpace, Filter blocks.
A unique feature of these blocks is that it is very convenient to quickly switch
between a continuous and a discrete block representation. Also, templates are provided
to quickly built-up standard controller structures.
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica_LinearSystems2.UsersGuide.GettingStarted\">Getting started</a>
     provides an overview of the Library in
     the <a href=\"modelica://Modelica_LinearSystems2.UsersGuide\">User's Guide</a>.</li>
<li><a href=\"modelica://Modelica_LinearSystems2.UsersGuide.ReleaseNotes\">Release Notes</a>
     summarizes the changes of new versions of this package.</li>
<li> <a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Contact\">Contact</a>
     gives the contact information for this library.</li>
<li> In an interactive environment, it is useful to run first the script
     &quot;_abbreviations.mos&quot; in directory
     &quot;Modelica_LinearSystems2/Resources/Scripts&quot;
     in order to set useful abbreviations: ss, tf, zp, poly, Complex, Plot, s, p, j.</li>
</ul>

<p>
It is planned to include this library in a future version of the
Modelica Standard Library.
</p>


<p>
<strong>Licensed by DLR under the 3-Clause BSD License</strong><br>
</p>

<table border=\"0\" cellpadding=\"2\" cellspacing=\"2\">
<tr>
  <td colspan=\"1\" rowspan=\"2\" style=\"vertical-align: middle;\">
    <img src=\"modelica://Modelica_LinearSystems2/Resources/Images/dlr_logo.png\">
  </td>
  <td style=\"vertical-align: top;\"><b>Copyright &copy; 2005-2012, DLR Institute of Robotics and Mechatronics</b></td>
</tr>
<tr>
  <td style=\"vertical-align: top;\"><b>Copyright &copy; 2012-2021, DLR Institute of System Dynamics and Control</b></td>
</tr>
</table>

<p>
<em>
This Modelica package is <u>free</u> software and
the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the
3-Clause BSD license, see the license conditions (including the
disclaimer of warranty) in the
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.The3clauseBSDLicense\">User's Guide</a>.
</em>
</p>

<p>
<strong>Modelica&reg;</strong> is a registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica_LinearSystems2;
