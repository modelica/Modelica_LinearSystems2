within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model Motor "A basic model of an electrical dc motor."

  Modelica.Electrical.Analog.Sources.SignalVoltage Vs annotation (Placement(
        transformation(
        origin={-70,0},
        extent={{10,-10},{-10,10}},
        rotation=90)));
  Modelica.Electrical.Analog.Basic.Ground G annotation (Placement(
        transformation(extent={{-80,-60},{-60,-40}}, rotation=0)));
  Modelica.Electrical.Analog.Basic.Resistor Ra(R=0.5) annotation (Placement(
        transformation(extent={{-60,30},{-40,50}}, rotation=0)));
  Modelica.Electrical.Analog.Basic.Inductor La(L=0.05) annotation (Placement(
        transformation(extent={{-20,30},{0,50}}, rotation=0)));
  Modelica.Electrical.Analog.Basic.EMF emf annotation (Placement(transformation(
          extent={{0,-10},{20,10}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput i_ref annotation (Placement(
        transformation(extent={{-108,-10},{-90,10}}, rotation=0)));
  Modelica.Mechanics.Rotational.Components.Inertia Jm(J=0.001) annotation (
      Placement(transformation(extent={{40,-10},{60,10}}, rotation=0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b annotation (
      Placement(transformation(extent={{90,-10},{110,10}}, rotation=0)));
equation
  connect(Ra.n, La.p) annotation (Line(points={{-40,40},{-20,40}}));
  connect(La.n, emf.p) annotation (Line(points={{0,40},{10,40},{10,10}}));
  connect(emf.flange, Jm.flange_a) annotation (Line(points={{20,0},{40,0}}));
  connect(Ra.p, Vs.p) annotation (Line(points={{-60,40},{-70,40},{-70,10}}));
  connect(Vs.n, emf.n)
    annotation (Line(points={{-70,-10},{-70,-20},{10,-20},{10,-10}}));
  connect(G.p, Vs.n) annotation (Line(points={{-70,-40},{-70,-10}}));
  connect(Jm.flange_b, flange_b) annotation (Line(points={{60,0},{100,0}}));
  connect(i_ref, Vs.v) annotation (Line(points={{-99,0},{-77,4.28626e-016}}));
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(
          extent={{-60,40},{60,-40}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-80,-80},{80,-100}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-38,-20},{-60,-80},{60,-80},{40,-20},{-40,-20},{-38,-20}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,0},{-60,0}}),
        Rectangle(
          extent={{60,8},{90,-8}},
          lineColor={160,160,164},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid),
        Text(extent={{-80,100},{80,60}}, textString="%name")}),
    Documentation(info="A basic model of an electrical dc motor.
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics));
end Motor;
