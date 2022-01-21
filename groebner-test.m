//filename:="6T15-[5,4,3]-51-42-33-g1.m";

import "reducecurve.m": SmallFunctions, PlaneModel, PlaneModelGroebner;

/*
K<nu> := ext<K|Polynomial(K, [-10, 0, 1])> where K is RationalField();
X := EllipticCurve([ 0, 0, 0, 1/46875*(-15884032*nu - 50229689), 1/52734375*(515007865088*nu + 1628597868634) ]);
KX<x,y> := FunctionField(X);
phi := (1/390625*(890044416*nu + 2814566400)*x^2 + 1/48828125*(4429534593024*nu + 14007418306560)*x + 1/6103515625*(6774503737196544*nu + 21422861828382720))/(x^6 + 1/125*(1216*nu + 3734)*x^5 + 1/9375*(-1047488*nu - 3310087)*x^4 + 1/10546875*(-116461946752*nu - 368285209652)*x^3 + 1/1318359375*(94825569207424*nu + 299864777724977)*x^2 + 1/2471923828125*(9827955791159757760*nu + 31078725043468227254)*x + 1/2780914306640625*(-117425818468119327543104*nu - 371333042468721457725271))*y + (1/390625*(-890044416*nu - 2814566400)*x^3 + 1/48828125*(-1761334198272*nu - 5569827840000)*x^2 + 1/30517578125*(-424143453578526720*nu - 1341259367958061056)*x + 1/3814697265625*(-1212572121351820476416*nu - 3834489730693837348864))/(x^6 + 1/125*(1216*nu + 3734)*x^5 + 1/9375*(-1047488*nu - 3310087)*x^4 + 1/10546875*(-116461946752*nu - 368285209652)*x^3 + 1/1318359375*(94825569207424*nu + 299864777724977)*x^2 + 1/2471923828125*(9827955791159757760*nu + 31078725043468227254)*x + 1/2780914306640625*(-117425818468119327543104*nu - 371333042468721457725271));

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
xs := SmallFunctions(PsQsRs, 2);
x_op := xs[4];
*/

F<g> := NumberField(Polynomial([8,-10,9,1,1]));
X := Curve(ProjectiveSpace(PolynomialRing(F, 2)));
KX1<x> := FunctionField(X);

tau := 2^38*3^17/23^3*(47323*g^3-1084897*g^2+7751*g-711002);
//_<x> := PolynomialRing(F);

P2 := (8*g^3+16*g^2-20*g+20)*x^2-(7*g^3+17*g^2-7*g+76)*x -
13*g^3+25*g^2-107*g+596;
P3 := 8*(31*g^3+405*g^2-459*g+333)*x^3+(941*g^3+1303*g^2-1853*g+1772)*x+85*g^3-
385*g^2+395*g-220;
P4 := 32*(4*g^3-69*g^2+74*g-49)*x^4+32*(21*g^3+53*g^2-68*g+58)*x^3-8*(97*g^3+95
*g^2-145*g+148)*x^2+8*(41*g^3-89*g^2-g+140)*x-123*g^3+391*g^2-93*g+3228;
P := P2^2*P3*P4^4/tau;

phi:=P;
Support(Divisor(phi));

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

xs := SmallFunctions(PsQsRs, 3);
for x_op in xs do
  pts, mults := Support(Divisor(x_op));
  printf "x_op has support\n%o,\n %o\n", pts, mults;
  time F_res := PlaneModel(phi, x_op);
  print "-------------------------------------";
end for;
//time F_gro := PlaneModelGroebner(phi, x_op);

/*
print "Computing image of map to AA2";
A2 := AffineSpace(F,2);
mp := map< X -> A2 | [x_op, phi]>;
time C := Image(mp);
F_img := DefiningEquation(C);
*/
