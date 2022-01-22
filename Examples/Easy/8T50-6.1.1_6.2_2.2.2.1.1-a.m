load "reducecurve.m";

// Belyi maps downloaded from the LMFDB on 21 January 2022.
// Magma code for Belyi map with label 8T50-6.1.1_6.2_2.2.2.1.1-a

// Geometric data

// Define the base field
R<T> := PolynomialRing(Rationals());
K<nu> := NumberField(R![-16, -6, 0, 1]);

// Define the curve
X := Curve(ProjectiveSpace(PolynomialRing(K, 2)));
// Define the map
KX<x> := FunctionField(X);
phi := (1/481555322138015200*(493227617344194*nu^2-8409605121382299*nu+625898368680545432)*x^8+1/15048603816812975*(-1009512635399912*nu^2+7989285079355496*nu-22347544826579504)*x^7+1/15048603816812975*(1651652380026974*nu^2-5893567739523146*nu+2371766838521296)*x^6)/(x^8+1/33042264685*(2050741150*nu^2-368305452*nu-25153405664)*x^7+1/15048603816812975*(-7232155966866374*nu^2+16837598850971684*nu+21397984515441104)*x^6+1/15048603816812975*(8444008412437728*nu^2-16131324482874912*nu-38060884899723776)*x^5+1/3009720763362595*(533908216336160*nu^2+1336939434642960*nu-10194320885419648)*x^4+1/15048603816812975*(-1107334029289984*nu^2-20478010887372800*nu+79470436567023616)*x^3+1/15048603816812975*(-619636611996672*nu^2+7969658983331328*nu-19530276227084288)*x^2+1/15048603816812975*(-5211564290170880*nu^2+9251408843759616*nu+26095240203796480)*x+1/15048603816812975*(2572636200759296*nu^2-2699836413274112*nu-19033152796737536));

print "phi has divisor";
Support(Divisor(phi));

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
printf "ramification points = %o\n", PsQsRs;

print "computing small functions supported at points above";
xs := SmallFunctions(PsQsRs, 2);
for x_op in xs do
  pts, mults := Support(Divisor(x_op));
  printf "x_op has support\n%o,\n %o\n", pts, mults;
  time F_res := PlaneModel(phi, x_op);
  print F_res;
  print Monomials(F_res);
  print "-------------------------------------";
end for;
