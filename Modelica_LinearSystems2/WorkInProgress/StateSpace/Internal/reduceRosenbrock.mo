within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
encapsulated function reduceRosenbrock
  "Algorithm to compress the generalized system matrix [A, B; C, D] to calculate the invariant zeros of a system"

  import Modelica;
  import Modelica.Math.Matrices.flipLeftRight;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  import Modelica.Utilities.Streams.print;

  input Real A[:,:];
  input Real B[:,:];
  input Real C[:,:];
  input Real D[:,:];

  output Real Ar[:,:];
  output Real Br[:,:];
  output Real Cr[:,:];
  output Real Dr[:,:];
  output Integer n;
  output Integer m;
  output Integer p;

protected
  Real A2[:,:];
  Real B2[:,:];
  Real C2[:,:];
  Real CC[:,:];
  Real Co[:,:];
  Real Cu[:,:];
  Real D2[:,:];
  Real DD[:,:];
  Real Mr[:,:];
  Real Vf[:,:];
  Real V[:,:];
  Real V2[:,:];
  Real R[:,:];
  Real tau[:];

  Integer nx=size(A, 1);
  Integer nu=size(B, 2);
  Integer ny=size(C, 1);

  Integer nue;
  Integer delta;
  Integer rho;
  Integer mue;
  Integer sigma;
  Integer j;
  Boolean stop;
  Boolean stop1 "reduction finished";
  Boolean stop2 "system has no zeros";
  Integer rankR;
  Real normA=Modelica.Math.Matrices.norm(A=A, p=1);
  Real eps=normA*1e-12;

algorithm
  if nx > 0 then

    A2 := A;
    B2 := B;
    C2 := C;
    D2 := D;
    stop := false;
    stop1 := false;
    stop2 := false;
    nue := nx;
    delta := 0;
    mue := ny;
    sigma := ny;
    j := 1;

    while not stop loop
      (V,R,tau,V2) := Matrices.QR( D2);

      rankR := 0;
//  !!!! rank has to be determined. In the case of ill conditioned systems svd should be used
      for i in 1:min(size(R, 1), size(R, 2)) loop
        if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
          rankR := rankR + 1;
        end if;
      end for;

//rankR:=Modelica.Math.Matrices.rank(R);

      DD := R[1:rankR, :];

      CC := transpose(V2)*C2;

      sigma := rankR;
      stop1 := size(CC,1) == rankR;

      if not stop1 then
        Cu := CC[sigma + 1:end, :];
        Co := CC[1:sigma, :];

        (V,R,tau,V2) := Matrices.QR( flipLeftRight(transpose(Cu)));
         Vf:=flipLeftRight(V2);

         rankR := 0;
//  !!!! rank determination
        for i in 1:min(size(R, 1), size(R, 2)) loop
          if abs(R[i, size(R, 2) - min(size(R, 1), size(R, 2)) + i]) > eps then
            rankR := rankR + 1;
          end if;
        end for;
//rankR:=Modelica.Math.Matrices.rank(R);

        rho := rankR;
        stop1 := rho == 0;
        stop2 := size(Cu,2) == rankR;

        if not stop1 and not stop2 then
          nue := size(Cu, 2) - rankR;
          mue := rho + sigma;
          delta := delta + rho;

          if sigma == 0 then
            Mr := [transpose(Vf)*A2*Vf,transpose(Vf)*B2];
          else
            Mr := [transpose(Vf)*A2*Vf,transpose(Vf)*B2; Co*Vf,DD];
          end if;

          A2 := Mr[1:nue, 1:nue];
          B2 := Mr[1:nue, nue + rho + 1:nue + rho + nu];
          C2 := Mr[nue + 1:nue + mue,1:nue];
          D2 := Mr[nue + 1:nue + mue,nue + rho + 1:nue + rho + nu];

          j := j + 1;
        end if;
       //not stop1 or not stop2

      end if;
       //if not stop1

      stop := stop1 or stop2 or j>3*nx;

    end while;

    if stop1 then
      Ar := A2;
      Br := B2;
      Cr := C2;
      Dr := D2;
      n := nue;
      p := sigma;
      m := nu;
    else
      n := 0;
      p := 0;
      m := 0;
      Ar := fill(0,0,0);
      Br := fill(0,0,0);
      Cr := fill(0,0,0);
      Dr := fill(0,0,0);
    end if;

  else
    n := 0;
    p := 0;
    m := 0;
    A2 := A;
    B2 := B;
    C2 := C;
    D2 := D;
  end if;

  annotation (Documentation(info="<html></html>"));
end reduceRosenbrock;
