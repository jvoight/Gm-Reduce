load "reducecurve.m";
Attach("monomials.m");

// Belyi maps downloaded from the LMFDB on 21 January 2022.
// Magma code for Belyi map with label 8T50-6.1.1_6.2_2.2.2.1.1-a

// Geometric data

// Define the base field
R<T> := PolynomialRing(Rationals());
K<nu> := NumberField(R![-16, -6, 0, 1]);
OK := Integers(K);

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

P := ideal< OK | OK![9973, 0, 0], OK![194, 1, 0] >;
printf "reducing mod P = %o\n", P;

X_FF, phi_FF := ReduceBelyiMap(X, phi, P);
RsandPs_FF := Support(Divisor(phi_FF));
RsandQs_FF := Support(Divisor(phi_FF-1));
PsQsRs_FF := SetToSequence(SequenceToSet(RsandPs_FF cat RsandQs_FF));
printf "ramification points = %o\n", PsQsRs_FF;

print "computing small functions supported at points above";
xs := SmallFunctions(PsQsRs, 2);
xs_FF := [];
for el in xs do
  _, x_FF := ReduceBelyiMap(X, el, P);
  Append(~xs_FF, x_FF);
end for;
for i := 1 to #xs do
  x_op := xs[i];
  x_op_FF := xs_FF[i];
  pts, mults := Support(Divisor(x_op));
  //printf "x_op has support\n%o,\n %o\n", pts, mults;
  print "computing model over finite field";
  F_res_FF := PlaneModel(phi_FF, x_op_FF);
  mons_FF := Monomials(F_res_FF);
  printf "%o monomials, max degree = %o\n", #mons_FF, Max([Degree(el) : el in mons_FF]);
  print "now computing in char 0";
  t0 := Cputime();
  F_res := PlaneModel(phi, x_op);
  t1 := Cputime();
  printf "computing plane model took %o\n", t1-t0;
  //print F_res;
  print "-------------------------------------";
end for;

//print "computing small functions supported at points above";
//xs := SmallFunctions(PsQsRs, 2);
//for x_op in xs do
//  pts, mults := Support(Divisor(x_op));
//  printf "x_op has support\n%o,\n %o\n", pts, mults;
//  print "computing plane model";
//  time F_res := PlaneModel(phi_FF, x_op);
//  print F_res;
//  print Monomials(F_res);
//  //  print "p-adic reduction";
//  //  time F_red_p := reducemodel_padic(F_res);
//  //  print "unit reduction";
//  //  time F_red_unit := reducemodel_units(F_red_p);
//  //  print F_red_unit;
//  //  print Monomials(F_red_unit);
//  print "-------------------------------------";
//end for;
