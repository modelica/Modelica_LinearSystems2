within Modelica_LinearSystems2.WorkInProgress.Tests.Conversion;
function convert_zp_2
  import Modelica_LinearSystems2.ZerosAndPoles;

protected
  ZerosAndPoles zp[15];

algorithm
  zp[1] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[2] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={1.1},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[3] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1={1.1,3.3},
    n2=fill(0, 0, 2), d2=fill(0, 0, 2));
  zp[4] := ZerosAndPoles(
    k=2,
    n1=fill(0,0),
    d1=fill(0,0),
    n2=fill(0, 0, 2), d2=[0,1]);
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
    d1={1,2,1,2,2.6},
    n2=[1, 1], d2=fill(0, 0, 2));
  zp[14] := ZerosAndPoles(
    n1=fill(0,0),
    d1=fill(0,0),
    n2=[1,2;-1,2], d2=[2,3;0,2]);

  zp[15]:= zp[1];
  for i in 2:14 loop
    zp[15]:=zp[15]*zp[i];
  end for;

  for i in 1:15 loop
    conv_zp2ss(zp[i],i);
  end for;

end convert_zp_2;
