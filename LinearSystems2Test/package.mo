within ;
package LinearSystems2Test "Library to test components of package LinearSystems2"
  extends Modelica.Icons.Package;

  annotation (
    preferredView="info",
    version="3.0.0",
    versionDate="2024-05-22",
    dateModified = "2024-05-22 19:40:00Z",
    uses(
      Modelica(version="4.0.0"),
      Modelica_LinearSystems2(version="3.0.0")),
    Icon(graphics={Text(
          extent={{-100,60},{100,-60}},
          textColor={0,0,0},
          textString="LS2")}),
    Documentation(info="<html>
<p>
This library provides models and functions to test components of
<strong>package Modelica_LinearSystems2</strong>.
</p>

<p>
Further development of this library should be performed in the following
way:
</p>

<ul>
  <li>
    Functions that are added to this library to test functions of the
    Modelica_LinearSystems2 Library should be called in
    <a href=\"modelica://LinearSystems2Test.testAllFunctions\">LinearSystems2Test.testAllFunctions</a>.
    The idea is that all test functions are called, when calling
    \"testAllFunctions()\".
  </li>
  <li>
    Models that are added to this library should have the annotation
    (with an appropriate StopTime):
    <blockquote><pre>
<strong>annotation</strong>(experiment(StopTime=1.1));
    </pre></blockquote>
    This gives the tool vendors the possibility to automatically identify
    the models that shall be simulated and, e.g., that shall be used in an automatic
    regression test.
  </li>
</ul>

<p>
Copyright &copy; 2024, DLR Institute of System Dynamics and Control.
</p>
<p>
For <strong>license</strong>, see
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.The3clauseBSDLicense\">3-clause BSD License</a>.
</p>

<p>
<strong>Modelica&reg;</strong> is a&nbsp;registered trademark of the Modelica Association.
</p>
</html>"));
end LinearSystems2Test;
