within ;
package ObsoleteLinearSystems2
  "Library that contains components from Modelica_LinearSystems2 Library 2.4.X that have been removed from version 3.0.0"

  package Math "Package of additional functions for Modelica.Math"
    extends Modelica.Icons.Package;

    operator record Complex "Record defining a Complex number"

      encapsulated function j "Obsolete imaginary unit function"
        import .Modelica;
        import .Complex;

        output Complex c = Modelica.ComplexMath.j "= sqrt(-1)";

        annotation(
          obsolete = "Obsolete function. See documentation for how to migrate imaginary unit j.",
          Inline=true,
          Documentation(info="<html>
<p>
A&nbsp;definition of complex number <code>c1</code> such as
</p>
<blockquote><pre>
Complex j = Modelica_LinearSystems2.Math.Complex.j();
Complex c1=2+3*j;
</pre></blockquote>
<p>
using Modelica_LinearSystems2 <strong>2.4.X</strong> shall be
automatically migrated now to (just the first line is affected)
</p>
<blockquote><pre>
Complex j = ObsoleteLinearSystems2.Math.Complex.j();
Complex c1=2+3*j;
</pre></blockquote>

<p>
The automatic conversion to proper <code>Modelica.ComplexMath.j</code> is
not possible since this is a&nbsp;complex constant. Therefore, you shall replace
the first line manually to
</p>
<blockquote><pre>
Complex j = Modelica.ComplexMath.j;
Complex c1=2+3*j;  // this stays unchanged
</pre></blockquote>
<p>
if <code>j</code> is further needed in the class.
In many classes, <code>j</code> is used only for the abovementioned complex
number definion. Then, the import of <code>j</code> can be completely omitted and
only the second line can be modified as follows
</p>
<blockquote><pre>
Complex c1 = Complex(2, 3);
</pre></blockquote>
</html>"));
      end j;

    end Complex;
  end Math;

  annotation (
    preferredView="info",
    version="3.0.0",
    versionDate="2024-06-21",
    dateModified = "2024-03-26 10:00:00Z",
    revisionId="$Format:%h %ci$",
    uses(
      Modelica(version="4.0.0"),
      Complex(version="4.0.0")),
    Documentation(info="<html>
<p>
This package contains functions and blocks from the Modelica_LinearSystems2 Library
version 2.4.X that are no longer available in the version 3.0.0.
The conversion script for version 2.4.X changes references in existing
user models automatically to the functions and blocks of package
ObsoleteLinearSystems2. The user should <strong>manually</strong> replace all
references to ObsoleteLinearSystems2 in his/her models to the functions
and blocks that are recommended in the documentation of the respective model.
</p>

<p>
In most cases, this means that a&nbsp;model with the name
\"ObsoleteLinearSystems2.XY\" should be renamed to \"Modelica_LinearSystems2.YZ\"
(version 3.0.0) and manually adaptated afterwards.
This usually requires some changes at the place where
the class is used (besides the renaming of the underlying class).
</p>

<p>
The models in ObsoleteLinearSystems2 are either not according to the
Modelica Language version 3.6 and higher, or the model was changed to get
a&nbsp;better design.
In all cases, an automatic conversion to the new implementation
was not feasible, since too complicated.
See also 
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.ReleaseNotes.Version_3_0_0\">Modelica_LinearSystems2.UsersGuide.ReleaseNotes.Version_3_0_0</a>.
</p>

<p>
In order to easily detect obsolete models and blocks, all of them are specially
marked in the icon layer with a red box. Additionally, an annotation &quot;obsolete&quot;
is provided.
</p>

<p>
<strong>Copyright</strong> &copy; 2024, DLR Institute of System Dynamics and Control
</p>

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
</html>"));
end ObsoleteLinearSystems2;
