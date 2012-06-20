within Modelica_LinearSystems2.Controller.UsersGuide;
class BlockProperties "Definition of block properties"

    annotation (Documentation(info="<html>
<p>Block <b>SampleClock</b> defines options for all components of the Controller library that are on the same or on a lower level as the sampleClock component. In every block, the default defined in SampleClock can be changed. The following options are available:</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/SampleClockMenu.png\"/> </p>
<p>With parameter <b>initTtype</b> the initialization of the blocks is defined. Default is steady state initialization (both for continuous and for discrete blocks). With parameter <b>blockType</b> the default block type is defined (Continuous or Discrete). If the Discrete block type is selected, the discretization method is defined via parameter <b>methodType</b> and the basic (periodic) sample time via parameter <b>sampleTime</b>. In every block, the input and output signals are sampled by &QUOT;sampleTime&QUOT; or an integer multiple of &QUOT;sampleTime&QUOT;. </p>
<p>In the <b>Advanced Options</b> menu of every block the following parameters are present:</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Controllers/AdvancedOptionsMenu.png\"/> </p>
<p>Parameters <b>blockType</b>, <b>methodType</b> and <b>initType</b> can be used to overwrite the default option defined in <b>SampleClock</b>. The <b>sampleFactor</b> defines the integer multiple of the base sample time of SampleClock that is used for the sampling of the input and output signals of this block. </p>
</html>"));
end BlockProperties;
