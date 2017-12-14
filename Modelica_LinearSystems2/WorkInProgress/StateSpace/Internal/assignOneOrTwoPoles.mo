within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
encapsulated function assignOneOrTwoPoles
  "Algorithm to assign p (p = 1 or 2) eigenvalues"

  import Modelica;
  import Modelica_LinearSystems2;
  import Complex;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real F[:,size(F, 1)] "system matrix of order p=1 or p=2";
  input Real G[size(F, 1),:] "control input matrix p rows";
  input Complex gamma[size(F, 1)];
  input Real tolerance=Modelica.Constants.eps;
  output Real K[:,size(F, 1)] "feedback matrix p columns";

protected
  Real Gamma[:,:];
  Integer rankGs;
  Real Fs[size(F, 1),size(F, 2)];
  Real Gs[size(G, 1),size(G, 2)];
  Real Gst[:,:]=transpose(G);
  Real Ks[:,size(F, 1)];
  Real c;
  Real s;
  Real r;
  Real z;
  Integer p=size(F,1);
  Real sigmaG[:];

  Real sumLambda = if p==1 then gamma[1].re else gamma[1].re+gamma[2].re;
  Real prodLambda = if p==1 then gamma[1].re else gamma[1].re*gamma[2].re-gamma[1].im*gamma[2].im;
  Real prodJ;

  Real j11;
  Real j12;
  Real j12_new;
  Real c0;
  Real c1;
  Real c3;
  Real c4;
  Real f_min_j12;
  Real der_f_min_j12;
  Integer inewt;
  Boolean den_nz;
  Boolean go_on;

  Real V[size(G, 2),size(G, 2)];
  Real U[size(F, 1),size(F, 2)];

  Real u1[:];
  Real u2[:];
  Integer i;
  Complex system_ev[:];

algorithm
  assert(size(F, 1) >= size(gamma, 1),
    "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F is of size ["
     + String(size(F, 1)) + "," + String(size(F, 1)) + "] and " + String(
    size(F, 1)) + " demanded assigned poles are expected. However, " +
    String(size(gamma, 1)) + " poles are given");
//assert(not Modelica.Math.Matrices.isEqual(G,zeros(size(G,1),size(G,2)),tolerance),"A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix ");
  if size(gamma, 1) == 1 then
    assert(gamma[1].im == 0, "\n In function StateSpace.Internal.assignOneOrTwoPoles() matrix F has size [" + String(size(F, 1)) + "," + String(size(F, 1)) +
      "], therefore, the demanded assigned pole must be real. However, the imaginary part is "
       + String(gamma[1].im));
  elseif abs(gamma[1].im) > 0 or abs(gamma[2].im) > 0 then
    assert(gamma[1].re == gamma[2].re and gamma[1].im == -gamma[2].im, "\nThe assigned pole pair given in function StateSpace.Internal.assignOneOrTwoPoles() must be conjungated complex. However, the poles are\npole1 = " + String(gamma[1]) + "\npole2 = " + String(gamma[2]) + ". \nTry\npole1 = " + String(gamma[1]) + "\npole2 = " + String(Modelica.ComplexMath.conj(gamma[1])) + "\ninstead");
  end if;

  if not Modelica.Math.Matrices.isEqual(
      G,
      zeros(size(G, 1), size(G, 2)),
      tolerance) then
    if size(G, 2) == 1 then
      V := [1];
      if size(G, 1) == 1 then
        U := [1];
      else
         // Givens
        r := sqrt(G[1, 1]^2 + G[2, 1]^2);
        c := G[1, 1]/r;
        s := G[2, 1]/r;
        U := [c,s; -s,c];
      end if;
      Gs := U*G;

      rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;
    else
     // size(G, 2)>1

      if size(G, 1) == 1 then // U=I, compute V by just one Householder transformation
        U := [1];
        u1 := cat(1, Vectors.householderVector(Gst[:, 1],
                     cat(1, {1}, zeros(size(G, 2) - 1))));// Householder vector
        Gst := Modelica_LinearSystems2.Math.Matrices.householderReflexion(Gst, u1);

        V := identity(size(G, 2)) - 2*matrix(u1)*transpose(matrix(u1))/(u1*u1);
        Gs := transpose(Gst);
        rankGs := if abs(Gs[1, 1]) > tolerance then 1 else 0;

      else
// systems with p==2 and m>1 are transformed by svd
        (sigmaG,U,V) := Modelica.Math.Matrices.singularValues(G);
        rankGs := 0;
        i := size(sigmaG, 1);
        while i > 0 loop
          if sigmaG[i] > 1e-10 then
            rankGs := i;
            i := 0;
          end if;
          i := i - 1;
        end while;
        Gs := zeros(p, size(G, 2));
        for i in 1:rankGs loop
          Gs[i, i] := sigmaG[i];
        end for;

      end if;
      V := transpose(V);
    end if;

// check controllability
    assert(not Modelica.Math.Matrices.isEqual(
      Gs,
      zeros(size(Gs, 1), size(Gs, 2)),
      tolerance), "A subsystem in StateSpace.Internal.assignOneOrTwoPoles() is not controllable");

    Ks := fill(
      0,
      rankGs,
      size(F, 1));
    Fs := U*F*transpose(U);

    if size(F, 1) == 1 then
      Ks := matrix((Fs[1, 1] - gamma[1].re)/Gs[1, 1]);
    else
      if rankGs == size(F, 1) then
        z := (Fs[1,1]+Fs[2,2] - sumLambda)/(Gs[1,1]*Gs[1,1]+Gs[2,2]*Gs[2,2]);
        Ks[1,1] := Gs[1,1]*z;//Minimum ks11 according to Frobenius norm
        Ks[2,2] := Gs[2,2]*z;//Minimum ks22 according to Frobenius norm
        j11 := Fs[1,1]-Gs[1,1]*Ks[1,1];
        prodJ := j11*(sumLambda-j11)-prodLambda;
        j12 := sqrt(abs(prodJ));// initial value of j12 element of the desiresd closed loop system matrix

// Newton iteration to find minimum ks12, ks21 according to Frobenius norm
        c0 := -prodJ*prodJ;
        c1 := prodJ*Fs[2,1];
        c4 := Gs[2,2]*Gs[2,2]/Gs[1,1]/Gs[1,1];
        c3 := -c4*Fs[1,2];

        inewt := 1;
        den_nz := true;
        go_on := true;
        while inewt<10 and den_nz and go_on loop
          f_min_j12 := c0 +j12*(c1 + j12*j12*(c3 + j12*c4));// minimum function
          der_f_min_j12 := c1 + j12*j12*(3*c3 + j12*4*c4);// derivation
          den_nz := abs(der_f_min_j12)>100*Modelica.Constants.eps;
          j12_new := if den_nz then j12 - f_min_j12/der_f_min_j12 else j12;
          go_on := abs(j12-j12_new)>100*Modelica.Constants.eps;
          j12 := j12_new;
          inewt := inewt+1;
        end while;

        Ks[1,2] := (Fs[1,2] - j12)/Gs[1,1];
        Ks[2,1] := if abs(j12)<100*Modelica.Constants.eps then (Fs[2,1] - sign(j12)*prodJ/100*Modelica.Constants.eps)/Gs[2,2] else (Fs[2,1] - prodJ/j12)/Gs[2,2];

      else

        Ks[1, 1] := (gamma[1].re + gamma[2].re - Fs[1, 1] - Fs[2, 2])/Gs[1, 1];
        Ks[1, 2] := Ks[1, 1]*Fs[2, 2]/Fs[2, 1] + (Fs[1, 1]*Fs[2, 2] - Fs[1, 2]*
          Fs[2, 1] - (gamma[1].re*gamma[2].re - gamma[1].im*gamma[2].im))/Fs[2,1]/Gs[1, 1];
        Ks := -Ks;
      end if;
    end if;

    K := V[:, 1:size(Ks, 1)]*Ks*U;

  else
    if p == 1 then
      Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, " + String(F[1, 1]) + " remains and " + String(gamma[1].re) + " cannot be realized");
    else
      system_ev :=Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(F);
      Modelica.Utilities.Streams.print("\n A subsystem (F, G) in StateSpace.Internal.assignOneOrTwoPoles() is not controllable, since G is equal to zero matrix. Therefore, K is set to zero matrix and the eigenvalues are retained.\n
      That is, " + String(system_ev[1].re) + (if abs(system_ev[1].im) > 0 then " + " else
              " - ") + String(system_ev[1].im) + "j and " + String(system_ev[2].re)
         + (if abs(system_ev[2].im) > 0 then " + " else " - ") + String(
        system_ev[2].im) + "j remain and " + String(gamma[1].re) + (if abs(
        gamma[1].im) > 0 then (if gamma[1].im > 0 then " + " else " - " +
        String(gamma[1].im) + "j") else "" + " and ") + String(gamma[2].re) + (
        if abs(gamma[2].im) > 0 then (if gamma[2].im > 0 then " + " else " - " +
        String(gamma[2].im) + "j") else "") + " cannot be realized");
    end if;
    K := zeros(size(G, 2), size(F, 1));
  end if;

end assignOneOrTwoPoles;
