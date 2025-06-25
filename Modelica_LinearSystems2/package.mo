within ;
package Modelica_LinearSystems2 "Modelica_LinearSystems2 (version 3.0.1) - Analysis, Synthesis and Modeling of Continuous and Discrete Linear Systems"

  extends Modelica.Icons.Package;

  constant String DataDir = Modelica.Utilities.Files.loadResource(
    "modelica://Modelica_LinearSystems2/Resources/Data/")
    "Absolute path to directory containing utility files for this package";


annotation (
  preferredView="info",
  uses(
    Complex(version="4.1.0"),
    Modelica(version="4.1.0"),
    DymolaCommands(version="1.20")),
  version="3.0.1",
  versionDate="2025-06-06",
  dateModified = "2025-06-06 15:44:03Z",
  revisionId="$Format:%h %ci$",
  conversion(
    from(
      version={"2.0", "2.1", "2.2", "2.3", "2.3.1", "2.3.2", "2.3.2", "2.3.3", "2.3.4"},
      to="2.3.5",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.4.mos"),
    from(
      version="2.3.5",
      to="2.4.0",
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.3.5.mos"),
    from(
      version={"2.4.0", "2.4.1"},
      script="modelica://Modelica_LinearSystems2/Resources/Scripts/Conversion/ConvertLinearSystems2_from_2.4.0.mos")),
  Documentation(
    info="<html>
<p>
Library <strong>Modelica_LinearSystems2</strong> is a&nbsp;Modelica package
providing different representations of linear, time invariant differential and
difference equation systems. For example, record
<a href=\"modelica://Modelica_LinearSystems2.StateSpace\">StateSpace</a>
defines a&nbsp;linear time invariant differential
equation system in state space form:
</p>
<blockquote><pre>
der(<strong>x</strong>) = <strong>A</strong> * <strong>x</strong> + <strong>B</strong> * <strong>u</strong>
    <strong>y</strong>  = <strong>C</strong> * <strong>x</strong> + <strong>D</strong> * <strong>u</strong>
</pre></blockquote>
<p>
Operators are overloaded to work conveniently with these system descriptions in an
interactive environment, e.g. to multiply transfer functions or to operate on complex numbers.
About 230 functions are provided to operate
on these data structures, e.g. to compute eigenvalues, zeros, step responses,
to design pole-placement and LQG controllers, to plot step responses, frequency responses,
eigenvalues, to convert between different description forms, or to
generate a&nbsp;linear system description by linearization of a&nbsp;Modelica model.
</p>

<p>
Furthermore, in subpackage
<a href=\"modelica://Modelica_LinearSystems2.Controllers\">Controllers</a>
about 20 input/output blocks of linear systems are provided that are
based on the different representation forms, e.g. PID, StateSpace or Filter blocks.
A&nbsp;unique feature of these blocks is that it is very convenient to quickly switch
between a&nbsp;continuous and a&nbsp;discrete block representation. Also, templates are provided
to quickly built-up standard controller structures.
</p>

<p>
For an introduction, have especially a&nbsp;look at:
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
     a&nbsp;built-in complex number type within Dymola.
     See also 
     <a href=\"modelica://Modelica_LinearSystems2.UsersGuide.GettingStarted.ComplexNumbers\">ComplexNumbers</a>
     in the <a href=\"modelica://Modelica_LinearSystems2.UsersGuide\">User's Guide</a>.</li>
</ul>

<p>
It is planned to include this library in a&nbsp;future version of the
Modelica Standard Library.
</p>

<p>
For <strong>copyright</strong> and BSD 3-clause <strong>license</strong>, see
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.The3clauseBSDLicense\">Copyright and License agreement</a>.
</p>

<p>
<strong>Modelica&reg;</strong> is a&nbsp;registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica_LinearSystems2;
