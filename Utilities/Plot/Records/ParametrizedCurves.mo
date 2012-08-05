within Modelica_LinearSystems2.Utilities.Plot.Records;
record ParametrizedCurves
  "Properties of a parameterized curve (displayed in a diagram)"
  extends Modelica.Icons.Record;

  Real X[:,:]
    "X=X(s), with X[j,i], where s[i] is the s-value and j is the s-branch"  annotation(Dialog);
  Real Y[size(X,1), size(X,2)]
    "Y=Y(s), with Y[j,i], where s[i] is the s-value and j is the s-branch"  annotation(Dialog);
  Real s[size(X,2)] "s[i] is the s-parametrization value for X and Y"  annotation(Dialog);
  String xName="x(s)" "Name of the x-variable (shown in tool tip)"  annotation(Dialog);
  String yName="y(s)" "Name of the y-variable (shown in tool tip)"  annotation(Dialog);
  String sName="s" "Name of the s-variable (shown in tool tip)"  annotation(Dialog);
  String heading="" "Heading displayed above diagram" annotation(Dialog);
  String legends[:]=fill("",0) "Legends of the curves" annotation(Dialog);
  Real heightRatio = 0.8 "Height of diagram = heightRatio*diagramWidth" annotation(Dialog);
  Boolean grid=true "True, if grid is shown" annotation(Dialog,  choices(__Dymola_checkBox=true));

  Boolean labelWithS=false "True, if s-values shall be shown along the curve"
       annotation(Dialog,  choices(__Dymola_checkBox=true));
  Modelica_LinearSystems2.Utilities.Plot.Records.CurveProperties curveProperties[:]=
      fill(Modelica_LinearSystems2.Utilities.Plot.Records.CurveProperties(),0)
    "Properties of the curves X[j,:] (if none given, a default is used; if only one curve property given, it is used for all curves)"
                                                                                                       annotation(Dialog);

  /* group "Axes" (Axes properties) */
  String xLabel=" " "String displayed at horizontal axis" annotation(Dialog(group="Axes"));
  String yLabel=" " "String displayed at vertical axis" annotation(Dialog(group="Axes"));
  Boolean logX = false "True, if logarithmic scale of x-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
  Boolean logY = false "True, if logarithmic scale of y-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
  Boolean uniformScaling = false
    "True, if same vertical and horizontal axis increment"
      annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));

  /* group "Legend" (Legend properties) */
  Boolean legend = false "True, if legend is shown" annotation(Dialog(group="Legend"),choices(__Dymola_checkBox=true));
  Boolean legendFrame=false "True, if frame around legend"
        annotation(Dialog(group="Legend"),   choices(__Dymola_checkBox=true));
  Boolean legendHorizontal=true
    "True, if horizontal legend (provided it is meaningful)"
        annotation(Dialog(group="Legend"),choices(__Dymola_checkBox=true));
  Modelica_LinearSystems2.Utilities.Plot.Types.LegendLocation legendLocation=
      Modelica_LinearSystems2.Utilities.Plot.Types.LegendLocation.Above
    "Legend placement" annotation(Dialog(group="Legend"));

end ParametrizedCurves;
