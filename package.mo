within ;
package Modelica_LinearSystems2 "Modelica_LinearSystems2 - Analysis, Synthesis and Modeling of Continuous and Discrete Linear Systems"

constant String DataDir=classDirectory() + "Extras/Data/"
  "Absolute path to directory containing utilitiy files for this package";













annotation (
  preferredView="info",
  uses(Modelica(version="3.1")),
   version="2.1",
   versionBuild=2,
   versionDate="2010-01-15",
   dateModified = "2010-01-15 09:27:58Z",
   revisionID="$Id::                                       $",
  Documentation(info="<html>
<p>
Library <b>Modelica_LinearSystems2</b> is a Modelica package
providing different representations of linear, time invariant differential and
difference equation systems. For example, record
<a href=\"Modelica://Modelica_LinearSystems2.StateSpace\">StateSpace</a>
defines a linear time invariant differential
equation system in state space form:
</p>
<pre>    <b>der</b>(x) = A * x + B * u
        y  = C * x + D * u
</pre>
<p>
Operators are overloaded to work conveniently with these system descriptions in an
interactive environment, e.g., to multiply transfer functions or to operate on complex numbers.
About 180 functions are provided to operate
on these data structures, e.g., to compute eigen values, zeros, step responses,
to design pole-placement and LQG controllers, to plot step responses, frequency responses,
eigen values, to convert between different description forms, or to
generate a linear system description by linearization of a Modelica model.
</p>

<p>
Furthermore, in sublibrary
<a href=\"Modelica://Modelica_LinearSystems2.Controller\">Controller</a>
about 20 input/output blocks of linear systems are provided that are
based on the different representation forms, e.g., PID, StateSpace, Filter blocks.
A unique feature of these blocks is that it is very convenient to quickly switch
between a continuous and a discrete block representation. Also, templates are provide
to quickly built-up standard controller structures.
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"Modelica://Modelica_LinearSystems2.UsersGuide.GettingStarted\">Getting started</a>
     provides an overview of the Library
     inside the <a href=\"Modelica://Modelica.UsersGuide\">Users Guide</a>.</li>
<li><a href=\"Modelica://Modelica_LinearSystems2.UsersGuide.ReleaseNotes\">Release Notes</a>
    summarizes the changes of new versions of this package.</li>
<li> <a href=\"Modelica://Modelica_LinearSystems2.UsersGuide.Contact\">Contact</a>
     gives the contact information for this library.</li>
<li> In an interactive environment, it is useful to run first the script
     &quot;_abbreviations.mos&quot; in directory
     &quot;Modelica_LinearSystems2\\Extras\\Scripts&quot;
     in order to set useful abbreviations: ss, tf, zp, poly, Complex, Plot, s, p, j.</li>
</ul>

<p>
It is planned to include this library in a future version of the
Modelica Standard Library. Note, the library is <u>not</u> backwards compatible
to the previous beta version 0.95, called \"Modelica_LinearSystems\", which was shipped
with previous versions of Dymola. Since the differences are too large, no conversion
scripts are provided, but different library names are used.
</p>


<p>
<b>Licensed by DLR under the Modelica License 2</b><br>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/dlr_logo.png\"  width=60 >
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
 <b>      Copyright &copy; 2005-2009, DLR Institute of Robotics and Mechatronics</b>
</p>


<p>
<i>This Modelica package is <u>free</u> software and
the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the
Modelica license 2, see the license conditions (including the
disclaimer of warranty)
<a href=\"Modelica://Modelica_LinearSystems2.UsersGuide.ModelicaLicense2\">here</a></u>
or at
<a href=\"http://www.Modelica.org/licenses/ModelicaLicense2\">
http://www.Modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</html>"));
end Modelica_LinearSystems2;
