within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model CraneWithEquations
  import SI = Modelica.SIunits;
  parameter SI.Mass m_crab=1000 "mass of crab";
  parameter SI.Mass m_load=4000 "mass of load";
  parameter SI.Length l=10 "length of rope";
  parameter SI.Acceleration g = 9.81 "Gravity acceleration";

  SI.Angle phi(start=0.5)
    "angle of rope with respect to gravity acceleration vector";
  SI.AngularVelocity w "= der(phi)";
  SI.Position s1;
  SI.Velocity v1 "horizontal velocity of crane";

  Modelica.Blocks.Interfaces.RealInput force
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y1 "Horziontal position of crab"
    annotation (Placement(transformation(extent={{100,50},{120,70}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y2 "horizontal position of load"
    annotation (Placement(transformation(extent={{100,-70},{120,-50}},
          rotation=0)));
equation
  y1=s1+l*sin(phi);
  y2=phi;

  der(phi)=w;
  der(w)=(-g*sin(phi)*(m_crab+m_load)-m_load*l*w^2*sin(phi)*cos(phi)-force*cos(phi))/(l*(m_crab+m_load*sin(phi)));
  der(s1)=v1;
  der(v1)=(m_load*g*cos(phi)*sin(phi)+m_load*l*w^2*sin(phi)+force)/(m_crab+m_load*sin(phi));

  annotation (uses(Modelica(version="3.1")),
    Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-80,70},{-40,40}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,40},{20,40}}, color={0,0,0}),
        Line(points={{-60,40},{0,-60}}, color={0,0,0}),
        Ellipse(
          extent={{-20,-40},{20,-80}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{96,48},{62,74}},
          lineColor={0,0,0},
          textString="y1"),
        Text(
          extent={{98,-74},{64,-48}},
          lineColor={0,0,0},
          textString="y2")}),
    experiment(
      StopTime=10,
      Interval=0.001,
      fixedstepsize=0.001,
      Algorithm="Euler"),
    experimentSetupOutput);
end CraneWithEquations;
