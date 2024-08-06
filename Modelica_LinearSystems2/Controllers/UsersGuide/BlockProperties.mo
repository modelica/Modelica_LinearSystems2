within Modelica_LinearSystems2.Controllers.UsersGuide;
class BlockProperties "Definition of block properties"
  extends Modelica.Icons.Information;

  annotation (
    Documentation(
      info="<html>
<p>
Block <a href=\"modelica://Modelica_LinearSystems2.Controllers.SampleClock\">SampleClock</a>
defines options for all components of the
<a href=\"modelica://Modelica_LinearSystems2.Controllers\">Controllers</a> library that are
on the same or on a&nbsp;lower level as the sampleClock component.
In every block, the default defined in SampleClock can be changed. The following options
are available:
</p>
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/SampleClockMenu.png\"/>
</div>
<p>
With parameter <code>initType</code> the initialization of the blocks is defined. Default
is steady state initialization (both for continuous and for discrete blocks). With parameter
<code>blockType</code> the default block type is defined (continuous or discrete).
If the discrete block type is selected, the discretization method is defined via parameter
<code>methodType</code> and the basic (periodic) sample time via parameter <code>sampleTime</code>.
In every block, the input and output signals are sampled by <code>sampleTime</code> or an
integer multiple of <code>sampleTime</code>.
</p>
<p>
In the <strong>Advanced Options</strong> menu of every block the following parameters are present:
</p>
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/AdvancedOptionsMenu.png\"/>
</div>
<p>
Parameters <code>blockType</code>, <code>methodType</code> and <code>initType</code> can be
used to overwrite the default option defined in <code>SampleClock</code>.
The <code>sampleFactor</code> defines the integer multiple of the base sample time of
SampleClock that is used for the sampling of the input and output signals of this block.
</p>
</html>"));
end BlockProperties;
