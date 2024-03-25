within Modelica_LinearSystems2.WorkInProgress.Controller;
block LimIntegratorReset
  "Integrator with limited values of the outputs, has external reset trigger"

  import Modelica.Blocks.Types.Init;
  extends Modelica.Blocks.Interfaces.SISO(y(start=y_start));

  Modelica.Blocks.Interfaces.RealInput lowerLimit
    "Connector of Real input signal"
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,-80})));
  Modelica.Blocks.Interfaces.RealInput upperLimit
    "Connector of Real input signal"
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,80})));

  parameter Real k(unit="1") = 1 "Integrator gain";
  parameter Modelica.Blocks.Types.Init initType=Modelica.Blocks.Types.Init.InitialState
    "Type of initialization (1: no init, 2: steady state, 3/4: initial output)"
    annotation(Evaluate=true, Dialog(group="Initialization"));
  parameter Boolean limitsAtInit=true
    "= false, if limits are ignored during initialization (i.e., der(y)=k*u)"
    annotation(Evaluate=true, Dialog(group="Initialization"));
  parameter Real y_start=0
    "Initial or guess value of output (must be in the limits lowerLimit .. upperLimit)"
    annotation (Dialog(group="Initialization"));

protected
  Boolean reset;
  Real limit;

initial equation
  if initType == Init.SteadyState then
    der(y) = 0;
  elseif initType == Init.InitialState or initType == Init.InitialOutput then
    y = y_start;
  end if;
equation
  if initial() and not limitsAtInit then
    der(y) = k*u;
    assert(y >= lowerLimit - 0.01*abs(lowerLimit) and y <= upperLimit + 0.01*
      abs(upperLimit), "LimIntegrator: During initialization the limits have been ignored.\n"
       + "However, the result is that the output y is not within the required limits:\n"
       + "  y = " + String(y) + ", lowerLimit = " + String(lowerLimit) + ", upperLimit = "
       + String(upperLimit));
  else
    der(y) = if y < lowerLimit or y > upperLimit then 0
       else k*u;
  end if;

  reset = y < lowerLimit or y > upperLimit;
  limit = if y < lowerLimit then lowerLimit else if y > upperLimit then upperLimit else y;
  when reset then
    reinit(y, limit);
  end when;

  annotation (
    Documentation(info="<html>
<p>This block extends <a href=\"modelica://Modelica.Blocks.Continuous.LimIntegrator\">Modelica.Blocks.Continuous.LimIntegrator</a> and adds an external reset trigger.  When the boolean input <strong>reset</strong> becomes true the integrator start is reset to <strong>y_start</strong>.</p>
</html>", revisions="<html>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={
        Rectangle(
          extent={{-100,-100},{100,100}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-74,-80},{26,20},{86,20}}, color={0,0,127}),
        Polygon(
          points={{-74,90},{-82,68},{-66,68},{-74,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{6,-10},{66,-70}},
          textColor={192,192,192},
          textString="I"),
        Polygon(
          points={{96,-80},{74,-72},{74,-88},{96,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-74,78},{-74,-90}}, color={192,192,192}),
        Line(points={{-84,-80},{88,-80}}, color={192,192,192})}),
    Diagram(graphics={
        Rectangle(extent={{-60,60},{60,-60}}),
        Line(points={{-100,0},{-60,0}}),
        Line(points={{60,0},{100,0}})}));
end LimIntegratorReset;
