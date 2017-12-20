within Modelica_LinearSystems2.Utilities.Types;
type StaircaseMethod = enumeration(
    QR "Apply staircase algorithm based on QR factorization",
    SVD "Apply staircase algorithm based on SVD")
  "Enumeration of methods for staircase algorithm"
    annotation (Evaluate=true,
    Icon(graphics={Ellipse(
        extent={{-100,100},{100,-100}},
        lineColor={255,0,128},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid), Text(
        extent={{-94,94},{94,-94}},
        lineColor={255,0,128},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid,
        textString="e")}));
