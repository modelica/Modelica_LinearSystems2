within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function drotg "Construct Givens plane rotation"
  extends Modelica.Icons.Function;

  input Real a;
  input Real b;
  output Real c;
  output Real s;
  output Real r=a;

external "Fortran 77" drotg(r,b,c,s) annotation(Library = {"lapack"});

  annotation (Documentation(info="Lapack documentation:

   Purpose
   =======

   construct givens plane rotation.

   =====================================================================  "));
end drotg;
