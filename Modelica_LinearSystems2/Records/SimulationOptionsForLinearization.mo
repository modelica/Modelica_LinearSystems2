within Modelica_LinearSystems2.Records;
record SimulationOptionsForLinearization
  "Options to define the simulation setup used for linearization"
  Boolean linearizeAtInitial=true
    "= true, if linearization at initial time; otherwise simulate until t_linearize"
     annotation (choices(checkBox=true));
  Modelica.Units.SI.Time t_start=0.0 "Start time of simulation"
    annotation (Dialog);
  Modelica.Units.SI.Time t_linearize=0.0
    "Simulate from t_start until t_linearize and then linearize, if linearizeAtInitial=false"
    annotation (Dialog(enable=not linearizeAtInitial));

  String method="Dassl" "Integration method, if linearizeAtInitial=false" annotation(Dialog(enable=not linearizeAtInitial),
      choices(
        choice="Lsodar" "Lsodar",
        choice="Dassl" "Dassl",
        choice="Euler" "Euler",
        choice="Rkfix2" "Rkfix2",
        choice="Rkfix3" "Rkfix3",
        choice="Rkfix4" "Rkfix4",
        choice="RadauIIa" "Radau IIa",
        choice="Esdirk23a" "Esdirk23a",
        choice="Esdirk34a" "Esdirk34a",
        choice="Esdirk45a" "Esdirk45a",
        choice="Dopri45" "Dopri45",
        choice="Dopri853" "Dopri853",
        choice="Sdirk34hw" "Sdirk34hw",
        choice="Cerk23" "Cerk23",
        choice="Cerk34" "Cerk34",
        choice="Cerk45" "Cerk45"));
  Real tolerance=1e-4 "Relative error tolerance, if linearizeAtInitial=false" annotation(Dialog(enable=not linearizeAtInitial));
  Real fixedStepSize=0.001
    "Step size for fixed step integrators, if linearizeAtInitial=false" annotation(Dialog(enable=not linearizeAtInitial));

  annotation (Icon(graphics={
        Rectangle(
          extent={{-100,20},{-72,-20}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-50,20},{-24,-20}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,20},{50,20},{50,50},{100,0},{50,-50},{50,-20},{0,-20},{0,20}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid)}));

end SimulationOptionsForLinearization;
