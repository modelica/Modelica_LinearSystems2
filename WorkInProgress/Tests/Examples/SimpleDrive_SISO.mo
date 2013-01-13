within Modelica_LinearSystems2.WorkInProgress.Tests.Examples;
model SimpleDrive_SISO
  parameter Modelica.SIunits.Radius r=0.5 "Radius of load";
  parameter Modelica.SIunits.Mass m=80 "Mass of load";
  DriveLib.Motor motor annotation (Placement(transformation(extent={{-10,-10},{
            10,10}}, rotation=0)));
  Modelica.Mechanics.Rotational.Components.IdealGear gearbox(ratio=100,
      useSupport=false) annotation (Placement(transformation(extent={{30,-10},{
            50,10}}, rotation=0)));
  Modelica.Mechanics.Rotational.Components.Inertia load(J=0.5*m*r*r)
    annotation (Placement(transformation(extent={{60,-10},{80,10}}, rotation=0)));
  Modelica.Mechanics.Rotational.Sensors.AngleSensor phiload annotation (
      Placement(transformation(
        origin={80,-30},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Interfaces.RealInput i_ref annotation (Placement(
        transformation(rotation=0, extent={{-110,10},{-90,30}})));
  Modelica.Blocks.Interfaces.RealOutput phi annotation (Placement(
        transformation(rotation=0, extent={{30,-110},{50,-90}})));
equation
  connect(gearbox.flange_b, load.flange_a)
    annotation (Line(points={{50,0},{60,0}}));
  connect(load.flange_b, phiload.flange)
    annotation (Line(points={{80,0},{80,-20}}));
  connect(motor.flange_b, gearbox.flange_a)
    annotation (Line(points={{10,0},{30,0}}));
  connect(i_ref, motor.inPort)
    annotation (Line(points={{-100,20},{-54.95,20},{-54.95,0},{-9.9,0}}));
  connect(phi, phiload.phi)
    annotation (Line(points={{40,-100},{40,-70.5},{80,-70.5},{80,-41}}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end SimpleDrive_SISO;
