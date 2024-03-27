within ;
package Modelica_LinearSystems2 "Modelica_LinearSystems2 (version 3.0.0-dev) - Analysis, Synthesis and Modeling of Continuous and Discrete Linear Systems"

  extends Modelica.Icons.Package;

  constant String DataDir = Modelica.Utilities.Files.loadResource(
    "modelica://Modelica_LinearSystems2/Resources/Data/")
    "Absolute path to directory containing utility files for this package";


annotation (
  preferredView="info",
  uses(
    Modelica(version="4.0.0"),
    DymolaCommands(version="1.17")),
  version="3.0.0-dev",
  versionDate="2024-06-21",
  dateModified = "2024-03-26 14:00:00Z",
  revisionId="$Format:%h %ci$",
  conversion(
    from(
      version={"2.0", "2.1", "2.2", "2.3", "2.3.1", "2.3.2", "2.3.2", "2.3.3", "2.3.4"},
      to="2.3.5",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.4.mos"),
    from(
      version="2.3.5",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.5.mos"),
    from(
      version={"2.4.0", "2.4.1"},
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.4.0.mos")),
  Documentation(
    info="<html>
<p>
Library <strong>Modelica_LinearSystems2</strong> is a Modelica package
providing different representations of linear, time invariant differential and
difference equation systems. For example, record
<a href=\"modelica://Modelica_LinearSystems2.StateSpace\">StateSpace</a>
defines a linear time invariant differential
equation system in state space form:
</p>
<blockquote><pre>
der(<strong>x</strong>) = <strong>A</strong> * <strong>x</strong> + <strong>B</strong> * <strong>u</strong>
    <strong>y</strong>  = <strong>C</strong> * <strong>x</strong> + <strong>D</strong> * <strong>u</strong>
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
<li> In an interactive environment, it is useful to run the script
     &quot;_abbreviations.mos&quot; in directory
     &quot;Modelica_LinearSystems2/Resources/Scripts&quot; first,
     in order to set useful abbreviations for e.g. ss, tf, zp, poly, j, etc.
     It is not necessary to import the package Complex since it is handled as
     a&nbsp;built-in complex number type within Dymola.</li>
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
  <td style=\"vertical-align: top;\"><strong>Copyright &copy; 2005-2012, DLR Institute of Robotics and Mechatronics</strong></td>
</tr>
<tr>
  <td style=\"vertical-align: top;\"><strong>Copyright &copy; 2012-2024, DLR Institute of System Dynamics and Control</strong></td>
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
