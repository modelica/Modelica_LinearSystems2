within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
encapsulated function reducedCtrSystem
  "calculate the controllable part of a SISO system"

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;
  import Modelica_LinearSystems2.Math.Complex;

  input StateSpace ss;

  output Modelica_LinearSystems2.Internal.StateSpaceR ssm1(
    redeclare Real A[size(ss.A, 1),size(ss.A, 2)],
    redeclare Real B[size(ss.B, 1),1],
    redeclare Real C[1,size(ss.C, 2)],
    redeclare Real D[size(ss.D, 1),size(ss.D, 2)])
    "controllable state space system";

protected
  Integer nx=size(ss.A, 1);
  Real Ah1[nx,size(ss.A, 2)];
  Real bh1[nx];
  Real ch1[nx];
  Real u[:] "householder vector";
  Integer ll;
  Integer r=1;
  Real maxa;
  Real normA=Modelica.Math.Matrices.norm(A=ss.A, p=1);
  Real eps = 2e-10;

algorithm
  if size(ss.C, 1) <> 1 or size(ss.B, 2) <> 1 then
    assert(size(ss.B, 2) == 1,
      "A SISO-system is expected as input\n but the number of inputs is "
       + String(size(ss.B, 2)) + " instead of 1");
    assert(size(ss.C, 1) == 1,
      " A SISO-system is expected as input\n but the number of outputs is "
       + String(size(ss.C, 1)) + " instead of 1");
  end if;

  Ah1 := ss.A;
  bh1 := ss.B[:, 1];
  ch1 := ss.C[1, :];

  if Modelica.Math.Vectors.length(bh1) > 0 then
    r := 1;
    if nx > 1 then

      u := Vectors.householderVector(bh1, cat(
            1,
            fill(0, nx - 1),
            {1}));

      Ah1 := Matrices.householderSimilarityTransformation(Ah1, u);
      bh1 := Vectors.householderReflexion_en(bh1, u);
      ch1 := Vectors.householderReflexion(ch1, u);
      bh1[1:nx - 1] := fill(0, nx - 1);

      ll := nx;
      maxa := max(abs(Ah1[1:ll - 1, ll]));

      while r <= nx - 1 and maxa > normA*eps loop
        u := cat(
              1,
              Vectors.householderVector(Ah1[1:ll - 1, ll], cat(
                1,
                fill(0, ll - 2),
                {1})),
              fill(0, nx - ll + 1));

        Ah1 := Matrices.Internal.hohoTrafoLowerHess(
              Ah1,
              u,
              r);

        ch1 := Vectors.householderReflexion(ch1, u);
        ll := ll - 1;
        maxa := max(abs(Ah1[1:ll - 1, ll]));

        r := r + 1;
      end while;

    end if;

    ssm1 := Modelica_LinearSystems2.Internal.StateSpaceR(
          A=Ah1,
          B=matrix(bh1),
          C=transpose(matrix(ch1)),
          D=ss.D,
          r=r);
    ssm1.r := r;
  end if;
end reducedCtrSystem;
