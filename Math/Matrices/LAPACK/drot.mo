within Modelica_LinearSystems2.Math.Matrices.LAPACK;
function drot "Applies a plane rotation"

  input Real x[:];
  input Real y[size(x,1)];
  input Real c;
  input Real s;
  input Integer incx=1;
  input Integer incy=1;

  output Real xr[size(x,1)]=x;
  output Real yr[size(x,1)]=y;

protected
  Integer n=size(x,1);

external "Fortran 77" drot(n,xr,incx,yr,incy,c,s) annotation(Library = {"lapack"});

  annotation (Documentation(info="
    Purpose
    =======
    applies a plane rotation.      "));
end drot;
