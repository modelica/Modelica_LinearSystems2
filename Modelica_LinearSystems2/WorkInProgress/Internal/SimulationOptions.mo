within Modelica_LinearSystems2.WorkInProgress.Internal;
record SimulationOptions "Options for simulation setup"
  Real startTime=0 "Start time of simulation" annotation (Dialog(
      group="Simulation Interval",
      __Dymola_label="Start time",
      __Dymola_absoluteWidth=15));
  Real stopTime=1 "Simulation stop time" annotation (Dialog(
      group="Simulation Interval",
      __Dymola_label="Stop time",
      __Dymola_absoluteWidth=15));
  Real outputInterval=0 "Distance between output points (if > 0)" annotation (
      Dialog(
      group="Output",
      __Dymola_label="Interval length",
      __Dymola_absoluteWidth=15));
  Integer numberOfIntervals=500
    "Number of intervals for output (if > 0 and Interval length == 0)"
    annotation (Dialog(
      group="Output",
      __Dymola_label="Number of intervals",
      __Dymola_absoluteWidth=15));
  String method "Integration method to be used" annotation(Dialog(
      group="Integration",
      __Dymola_label="Algorithm",
      __Dymola_absoluteWidth=15),
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
  Real tolerance=1e-3 "Relative error tolerance" annotation (
      Dialog(
      group="Integration",
      __Dymola_label="Tolerance",
      __Dymola_absoluteWidth=15));
  Real fixedStepSize=0 "Step size for fixed step integrators" annotation (
      Dialog(
      group="Integration",
      __Dymola_label="Fixed integrator step",
      __Dymola_absoluteWidth=15));
  annotation (__Dymola_Protection(
      hideFromBrowser=false,
      allowDuplicate=true,
      showDiagram=true,
      showText=true,
      showVariables=true,
      showDiagnostics=true,
      showStatistics=true), Icon(graphics={
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
          points={{0,20},{50,20},{50,50},{100,0},{50,-50},{50,-20},{0,-20},{0,
              20}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid)}));
end SimulationOptions;
