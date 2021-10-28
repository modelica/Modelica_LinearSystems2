within Modelica_LinearSystems2.Controller.UsersGuide;
class BlockProperties "Definition of block properties"

    annotation (Documentation(info="<html>
<p>Block <strong>SampleClock</strong> defines options for all components of the Controller library that are on the same or on a lower level as the sampleClock component. In every block, the default defined in SampleClock can be changed. The following options are available:</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/SampleClockMenu.png\"/> </p>
<p>With parameter <strong>initTtype</strong> the initialization of the blocks is defined. Default is steady state initialization (both for continuous and for discrete blocks). With parameter <strong>blockType</strong> the default block type is defined (Continuous or Discrete). If the Discrete block type is selected, the discretization method is defined via parameter <strong>methodType</strong> and the basic (periodic) sample time via parameter <strong>sampleTime</strong>. In every block, the input and output signals are sampled by &quot;sampleTime&quot; or an integer multiple of &quot;sampleTime&quot;. </p>
<p>In the <strong>Advanced Options</strong> menu of every block the following parameters are present:</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/AdvancedOptionsMenu.png\"/> </p>
<p>Parameters <strong>blockType</strong>, <strong>methodType</strong> and <strong>initType</strong> can be used to overwrite the default option defined in <strong>SampleClock</strong>. The <strong>sampleFactor</strong> defines the integer multiple of the base sample time of SampleClock that is used for the sampling of the input and output signals of this block. </p>
</html>"));
end BlockProperties;
