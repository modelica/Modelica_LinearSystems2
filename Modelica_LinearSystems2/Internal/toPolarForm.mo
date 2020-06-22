within Modelica_LinearSystems2.Internal;
encapsulated function toPolarForm
  "Transform a vector of complex values defined as a real Zeros matrix to polar form"
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Internal;

  input Real Zeros[:,2]
    "Zeros as Real matrix (first column: real, second column imaginary values)";
  output Real A[size(Zeros,1)] "Amplitudes of Zeros";
  output Modelica.Units.SI.Angle phi[size(Zeros, 1)] "Angles of zeros";
  output Boolean oneZeroAmplitude
    "=true: at least one element of Zeros has zero amplitude";
protected
  Real r_abs;
  Real i_abs;
algorithm
  oneZeroAmplitude :=false;
  for i in 1:size(Zeros,1) loop
    r_abs :=abs(Zeros[i,1]);
    i_abs :=abs(Zeros[i,2]);
    if r_abs == 0 and i_abs == 0 then
      A[i] :=0;
      phi[i] :=0;
      oneZeroAmplitude :=true;
    else
      A[i]   :=if r_abs > i_abs then r_abs*sqrt(1 + (i_abs/r_abs)^2) else i_abs*sqrt(1 + (r_abs/i_abs)^2);
      phi[i] :=atan2(Zeros[i,2], Zeros[i,1]);
    end if;
  end for;
end toPolarForm;
