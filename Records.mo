within Modelica_LinearSystems2;
package Records "Records, especially to build menus"
extends Modelica.Icons.Package;

  record ParameterVariation
    "Define variation of one parameter in a given range and optionally select the parameter from a translated model"
    String Name "Name of parameter" annotation (Dialog);
    Real Min "Minimum value of parameter" annotation (Dialog);
    Real Max "Maximum value of parameter" annotation (Dialog);
    Integer nVar(min=2) = 10 "Number of parameter variations (min=2)" annotation (Dialog);
    String Unit="" "Unit of parameter" annotation (Dialog);
    annotation (
    Dialog(__Dymola_importDsin(onlyStart=true,
      fields(Name=initialName,
             Min=initialValue.minimum,
             Max=initialValue.maximum,
             Unit=initialValue.unit))),
    Icon(graphics={
          Rectangle(
            extent={{-100,-30},{100,-90}},
            lineColor={0,0,0},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-100,-30},{-100,44}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{100,-32},{100,44}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-100,40},{100,40}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-60,68},{-100,40},{-60,12}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{60,68},{100,40},{60,12}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end ParameterVariation;

  record SimulationOptionsForLinearization
    "Options to define the simulation setup used for linearization"
    String method="Dassl" "Integration method" annotation(Dialog,
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
    Real tolerance=1e-4 "Relative error tolerance" annotation(Dialog);
    Real fixedStepSize=0.001 "Step size for fixed step integrators" annotation(Dialog);

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
end Records;
