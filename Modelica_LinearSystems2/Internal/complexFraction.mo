within Modelica_LinearSystems2.Internal;
function complexFraction
  "Return z = n[1]*n[2]*...n[end]/(d[1]*d[2]*....*d[end])"
  input Complex n[:];
  input Complex d[:];
  output Complex z "= n[1]*n[2]*...n[end]/(d[1]*d[2]*....*d[end])";
  output Integer info
    "= 0/1/2 success/infinity(z is a large value)/indefinite (0/0; z=0)";
protected
  Integer nn = size(n,1);
  Integer nd = size(d,1);
  Real n_abs[nn];
  Real n_abs2[nn];
  Real n_phi[nn];
  Real d_abs[nd];
  Real d_abs2[nd];
  Real d_phi[nd];
  Real wr;
  Real wi;
  Real wr_abs;
  Real wi_abs;
  Boolean n_zero = false;
  Boolean d_zero = false;
  Real z_abs;
  Real z_phi;
  constant Real pi = Modelica.Constants.pi;
  constant Real pi2 = 2*pi;
algorithm
  // Transform elements of n to polar form
  for i in 1:nn loop
    wr :=n[i].re;
    wi :=n[i].im;
    wr_abs :=abs(wr);
    wi_abs :=abs(wi);
    if wr_abs == 0 and wi_abs == 0 then
      n_abs[i] :=0;
      n_phi[i] :=0;
      n_zero :=true;
    else
      n_abs[i] :=if wr_abs > wi_abs then wr_abs*sqrt(1 + (wi/wr)^2) else wi_abs*sqrt(1 + (wr/wi)^2);
      n_phi[i] :=atan2(wi, wr);
    end if;
  end for;

  // Transform elements of d to polar form
  for i in 1:nd loop
    wr :=d[i].re;
    wi :=d[i].im;
    wr_abs :=abs(wr);
    wi_abs :=abs(wi);
    if wr_abs == 0 and wi_abs == 0 then
      d_abs[i] :=0;
      d_phi[i] :=0;
      d_zero :=true;
    else
      d_abs[i] :=if wr_abs > wi_abs then wr_abs*sqrt(1 + (wi/wr)^2) else wi_abs*sqrt(1 + (wr/wi)^2);
      d_phi[i] :=atan2(wi, wr);
    end if;
  end for;

  // Handle zeros in n
  if n_zero then
    z :=Complex(0, 0);
    info :=if d_zero then 2 else 0;
    return;
  end if;

  // Compute angle of fraction
  z_phi :=0;
  if nn <= nd then
    for i in 1:nn loop
      z_phi := z_phi + n_phi[i] - d_phi[i];
    end for;

    for i in nn+1:nd loop
      z_phi := z_phi - d_phi[i];
    end for;
  else
    for i in 1:nd loop
      z_phi := z_phi + n_phi[i] - d_phi[i];
    end for;

    for i in nd+1:nn loop
      z_phi := z_phi + n_phi[i];
    end for;
  end if;

  // Compute absolute value (avoid overflow)
  if d_zero then
    info :=1;
    z_abs :=Modelica.Constants.inf;
  else
    info :=0;
    if nn > 0 then
      n_abs2 :=Modelica.Math.Vectors.sort(n_abs, ascending=false);
    end if;
    if nd > 0 then
      d_abs2 :=Modelica.Math.Vectors.sort(d_abs, ascending=false);
    end if;
    z_abs := 1;

    if nn <= nd then
      for i in 1:nn loop
        z_abs :=z_abs*(n_abs2[i]/d_abs2[i]);
      end for;
      for i in nn+1:nd loop
        z_abs :=z_abs/d_abs2[i];
      end for;
    else
      for i in 1:nd loop
        z_abs :=z_abs*(n_abs2[i]/d_abs2[i]);
      end for;
      for i in nd+1:nn loop
        z_abs :=z_abs*n_abs2[i];
      end for;
    end if;
  end if;

  // Transform polar in Cartesian description
  z.re :=z_abs*cos(z_phi);
  z.im :=z_abs*sin(z_phi);
end complexFraction;
