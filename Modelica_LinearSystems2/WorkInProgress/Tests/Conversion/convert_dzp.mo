within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function convert_dzp
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;
protected
  ZerosAndPoles zp[15];
  Modelica_LinearSystems2.DiscreteZerosAndPoles dzp[
                            15];

algorithm
  zp[1] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[2] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={0},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[3] := ZerosAndPoles(
    k=10,
    n1=fill(0,0),
    d1={0,0},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[4] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=[0,0]);
  zp[5] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={1,2},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[6] := ZerosAndPoles(
    k=2,
    n1={-1},
    d1={1,2},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[7] := ZerosAndPoles(
    k=2,
    n1={-1,3},
    d1={1,2},
    n2=fill(0, 0, 2), d2=fill(0,0,2));
  zp[8] := ZerosAndPoles(
    k=2,
    n1={0,-1,2,-3,4},
    d1={1,3,2.5,3.5,1},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[9] := ZerosAndPoles(
    k=2,
    n1={1,2},
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=[2,2]);
  zp[10] := ZerosAndPoles(
    n1={1,2,1.5},
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=[2,2;1,3]);
  zp[11] := ZerosAndPoles(
    k=2,
    n1={1,2,1.5},
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=[2,2;1,3;-1,1.5]);
  zp[12] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={1,2},
    n2=[1, 1], d2=fill(0, 0, 2));
  zp[13] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={1,2,1,2,0},
    n2=[1, 1], d2=fill(0, 0, 2));
  zp[14] := ZerosAndPoles(
    n1=fill(0,0),
    d1=fill(0,0),
    n2=[1,2;-1,2], d2=[2,3;0,0]);

  for i in 1:14 loop
    dzp[i] :=DiscreteZerosAndPoles(
      zp[i],
      Ts=0.1,
      method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact);
  end for;

  dzp[15]:= dzp[1];
  for i in 1:10 loop
    dzp[15]:=dzp[15]*dzp[i];
  end for;

  for i in 1:15 loop
    conv_dzp2dss(dzp[i],i);
  end for;

end convert_dzp;
