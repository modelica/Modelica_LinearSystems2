within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function drotg "Construct Givens plane rotation"

  input Real a;
  input Real b;
  output Real c;
  output Real s;
  output Real r=a;

external "Fortran 77" drotg(r,b,c,s) annotation(Library = {"lapack"});

  annotation (Documentation(info="
    Purpose   
    =======   
    construct givens plane rotation.     "));
end drotg;
