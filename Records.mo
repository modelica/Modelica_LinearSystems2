within Modelica_LinearSystems2;
package Records "Records, especially to build menus"
extends Modelica.Icons.Package;

  record SetParameter "Set value of one parameter from a translated model"

    String Name "Name of parameter" annotation (Dialog);
    Real Value "Value of parameter" annotation (Dialog);

    annotation (
    Dialog(__Dymola_importDsin(button="select"
            "Select the model parameters to be included in this table",
                                              onlyStart=true,
      fields(Name=initialName,
             Value=initialValue.value))),
    Icon(graphics={
          Rectangle(
            extent={{-100,-30},{100,-90}},
            lineColor={0,0,0},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-40,-30},{-40,80}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{20,28},{-20,0},{20,-28}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={-40,-10},
            rotation=90)}));

  end SetParameter;

  record ParameterVariation
    "Define variation of one parameter in a given range and optionally select the parameter from a translated model"

    String Name "Name of parameter" annotation (Dialog);
    Modelica_LinearSystems2.Types.Grid grid = Modelica_LinearSystems2.Types.Grid.OneValue
      "Definition of parameter grid"    annotation (Dialog);

    Real Value=0 "Value of parameter" annotation (Dialog(group="if grid = OneValue"));
    Real Min=-1e100 "Minimum value of parameter"
                                   annotation(Dialog(group="if grid = Equidistant or Logarithmic"));
    Real Max=1e100 "Maximum value of parameter"
                                   annotation(Dialog(group="if grid = Equidistant or Logarithmic"));
    Integer nPoints(min=2) = 11
      "Number of parameter values in the range Min .. Max"    annotation(Dialog(group="if grid = Equidistant or Logarithmic"));

    annotation (
    Dialog(__Dymola_importDsin(button="select"
            "select the model parameters to be included in this table",                                   onlyStart=true,
      fields(Name=initialName,
             Value=initialValue.value,
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
    Boolean linearizeAtInitial=true
      "= true, if linearization at inital time; otherwise simulate until t_linearize"
       annotation (Dialog,choices(__Dymola_checkBox=true));
    Modelica.SIunits.Time t_start=0.0 "Start time of simulation" annotation(Dialog);
    Modelica.SIunits.Time t_linearize=0.0
      "Simulate from t_start until t_linearize and then linearize, if linearizeAtInitial=false"
                              annotation(Dialog(enable=not linearizeAtInitial));

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
    Real tolerance=1e-4 "Relative error tolerance, if linearizeAtInitial=false"
                                                                                annotation(Dialog(enable=not linearizeAtInitial));
    Real fixedStepSize=0.001
      "Step size for fixed step integrators, if linearizeAtInitial=false"                        annotation(Dialog(enable=not linearizeAtInitial));

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
