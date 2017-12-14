within Modelica_LinearSystems2.WorkInProgress.StateSpace.Internal;
function assignPolesSI_rq
  "RQ implementation of a recursiv single-input pole assignment problem"
  extends Modelica.Icons.Function;

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.Math.Matrices;
  import Complex;
  import matMul = Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.matMatMul;
  import Modelica_LinearSystems2.Math.ComplexAdvanced.Matrices.matVecMul;
  import Modelica_LinearSystems2.Math.ComplexAdvanced.Internal.C_transpose;
  import Re = Modelica.ComplexMath.real;
  import Im = Modelica.ComplexMath.imag;
  import Modelica.Utilities.Streams.print;

  input StateSpace ss;
  input Complex gamma[size(ss.A, 1)];

  output Real K[1,size(ss.A, 1)] "feedback matrix";
  output Real S[size(ss.A, 1),size(ss.A, 1)] "Closed loop System matrix";
  output Complex po[size(ss.A, 1)] "poles of the closed loop system";
protected
  Boolean isCntrbl "is set to true if the system is controllable";
  Modelica_LinearSystems2.Internal.StateSpaceR ssr "system in Hessenberg form";

  Integer nx=size(ss.A, 1);
  Complex H[size(ss.A, 1),size(ss.A, 1)];
  Complex Q[size(ss.A, 1),size(ss.A, 1)];
  Complex R[size(ss.A, 1),size(ss.A, 1)];
  Complex Qi[:,:];
  Complex Ri[:,:];
  Complex Hi[:,:];
  Complex Hlr[:,:];
  Complex Qn[size(ss.A, 1),size(ss.A, 1)];
  Complex rho;

  Real alpha=1.0;
  Real beta;
  Integer i;
  Integer ni;

  Real P[:,:];

  Complex ev[:];
  Complex RQ[:,:];
  Complex tau[:];

algorithm
  //reduction to controller-Hessenberg form
  (isCntrbl,ssr,P) := StateSpace.Internal.staircaseSVD(ss);

  assert(isCntrbl, "Poles cannot be placed since system is not controllable");

  H := Complex(1)*ssr.A;

  if nx == 1 then
    K[1, 1] := (ss.A[1, 1] - Re(gamma[1]))/ss.B[1,1];
  else
    beta := ssr.B[1,1];
    for i in 2:nx loop
      alpha := alpha*Re(H[i, i - 1]);
    end for;

    for i in 1:nx loop
      H[i, i] := H[i, i] - gamma[1];
    end for;

    (R,Q) := Modelica_LinearSystems2.WorkInProgress.Math.Matrices.C_RQ(
                           H);
    H := matMul(Q, R);

    for i in 1:nx loop
      H[i, i] := H[i, i] + gamma[1];
    end for;

    rho := R[nx, nx];
    Hi := H;

    for i in 2:nx loop
      ni := size(Hi, 1);
      for ii in 1:ni loop
        Hi[ii, ii] := Hi[ii, ii] - gamma[i];
      end for;

//       (Ri,Qi) :=  Matrices.C_RQ(Hi);
//       Hlr := Qi*Ri;
// replaced by :########
      (RQ,tau) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgerq2(
                                         Hi);
      Ri := fill(Complex(0),ni,ni);
      for ii in 1:ni loop
        for iii in ii:ni loop
          Ri[ii,iii] := RQ[ii,iii];
        end for;
      end for;
      Qi := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zungrq(
                                   RQ,tau);
      Hlr := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zunmrq(
                                    Ri,RQ,tau,true);
// ##########

      for ii in 1:ni loop
        Hlr[ii, ii] := Hlr[ii, ii] + gamma[i];
      end for;
      Hi := Hlr[2:ni, 2:ni];

      if i > 2 then
        Qn := fill(Complex(0), nx, nx);
        for ii in 1:i - 2 loop
          Qn[ii, ii] := Complex(1);
        end for;
        for ii in 1:nx-i+2 loop
          for iii in 1:nx-i+2 loop
            Qn[ii+i-2, iii+i-2] := Qi[ii, iii];
          end for;
        end for;
        Q := matMul(Qn,Q);
      else
        Q := matMul(Qi,Q);
      end if;

      rho := rho*Ri[nx + 2 - i, nx + 2 - i];
    end for;

    K[1, :] := Q[nx, :].re;
    K := (rho.re/alpha/beta)*K*P;
    end if;
    S := ss.A-ss.B*K;
    po := Modelica_LinearSystems2.Math.ComplexAdvanced.eigenValues(S);
  Modelica_LinearSystems2.Math.ComplexAdvanced.Vectors.print("ev", ev);
end assignPolesSI_rq;
