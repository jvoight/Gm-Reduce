/*
  AttachSpec("code/spec_database");
  load "../../Gm Reduce/Code/reducecurve.m";
  filename:="6T15-[5,4,3]-51-42-33-g1.m";
  X:=BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
  phi:=BelyiDBRead(filename)`BelyiDBBelyiMaps[1];
  X; phi;
*/

AttachSpec("spec");
K<nu> := ext<K|Polynomial(K, [-10, 0, 1])> where K is RationalField();
X := EllipticCurve([ 0, 0, 0, 1/46875*(-15884032*nu - 50229689), 1/52734375*(515007865088*nu + 1628597868634) ]);
KX<x,y> := FunctionField(X);
phi := (1/390625*(890044416*nu + 2814566400)*x^2 + 1/48828125*(4429534593024*nu + 14007418306560)*x + 1/6103515625*(6774503737196544*nu + 21422861828382720))/(x^6 + 1/125*(1216*nu + 3734)*x^5 + 1/9375*(-1047488*nu - 3310087)*x^4 + 1/10546875*(-116461946752*nu - 368285209652)*x^3 + 1/1318359375*(94825569207424*nu + 299864777724977)*x^2 + 1/2471923828125*(9827955791159757760*nu + 31078725043468227254)*x + 1/2780914306640625*(-117425818468119327543104*nu - 371333042468721457725271))*y + (1/390625*(-890044416*nu - 2814566400)*x^3 + 1/48828125*(-1761334198272*nu - 5569827840000)*x^2 + 1/30517578125*(-424143453578526720*nu - 1341259367958061056)*x + 1/3814697265625*(-1212572121351820476416*nu - 3834489730693837348864))/(x^6 + 1/125*(1216*nu + 3734)*x^5 + 1/9375*(-1047488*nu - 3310087)*x^4 + 1/10546875*(-116461946752*nu - 368285209652)*x^3 + 1/1318359375*(94825569207424*nu + 299864777724977)*x^2 + 1/2471923828125*(9827955791159757760*nu + 31078725043468227254)*x + 1/2780914306640625*(-117425818468119327543104*nu - 371333042468721457725271));

f:=UnreducedModel(phi);



prec := 0;
K := BaseRing(Parent(f));
if prec eq 0 then
  //wild guess imprecise
  prec:=Floor(Sqrt(Degree(K)))*100;
end if;

//u := Parent(fuv).1;
//v := Parent(fuv).2;
ZK := Integers(K);
r,s:=Signature(K);
k:=RealField(prec);
//Rx3<x1,x2,x3>:=PolynomialRing(Rationals(),3);

variables:=[ Parent(f).i : i in [1..#Names(Parent(f))] ];
var_size:=#variables;
ZK := Integers(K);

inf_places:=InfinitePlaces(K);
assert #inf_places eq r+s;
phi:=function(x);
  return [ Log(k!Abs(Evaluate(x,v : Precision:=prec))) : v in inf_places ];
end function;

mexps := [ Exponents(m) : m in Monomials(f) ];
coefs:=Coefficients(f);
//assert &+[ coefs[i]*(u^mexps[i,1])*v^mexps[i,2] : i in [1..#mexps] ] eq fuv;

UK,mUK:=UnitGroup(K);
k := RealField(prec);
//UU:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) ];
UU:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) and k!0 notin phi(K!(mUK(eps)))  ];
RRUU := [ phi(eps) : eps in UU ];
L0gens:= [ &cat[ u : i in [1..#mexps] ] : u in RRUU ];
Lambda:= Lattice((r+s)*#mexps, L0gens[1]);


constants := [];
abs_coef := [];

for n in [1..#mexps] do
  alpha_norm := Log(k!Abs(Norm(coefs[n])))/(r+s);
  log_coef:= phi(coefs[n]);
  tuple:=[];

  for m in [1..r+s] do
    if m le r then
      const:= Log(Abs(Evaluate(coefs[n], inf_places[m] : Precision:=prec))) - Log(k!Abs(Norm(coefs[n])))/(r+s);
    else
      const:= Log(Abs(Evaluate(coefs[n], inf_places[m] : Precision:=prec))) - Log(k!Abs(Norm(coefs[n])))/(2*(r+s));
    end if;
    Append(~tuple, const);

    etas:= [ phi(eps)[m] : eps in UU ];
    lhs:=&cat[ [ eta*a : a in mexps[n] ] cat [eta] : eta in etas ];
    Append(~abs_coef, lhs);

  end for;
  Append(~constants,tuple);
end for;

average:= [ &+[ (A[i])/#mexps : A in constants ] : i in [1..r+s] ];


V := VectorSpace(k,r+s);
V0:= VectorSpace(k,(r+s)*#mexps);



L2coef:= [ &+[I[i] : I in mexps] : i in [1..var_size] ] cat [#mexps];
L2pols<[X]> := PolynomialRing(k,3*(r+s));
quadform_pol:= &+[ (L2coef[1]*X[i] + L2coef[1]*X[r+s+i] + L2coef[1]*X[2*(r+s)+i])^2 : i in [1..r+s] ];
L2mat:=SymmetricMatrix(quadform_pol);

constants := &cat constants;
beta:=V0!constants;

//Lambda :=

  constants:=Matrix(k,constants);
  abs_coef:=Matrix(k,abs_coef);

  L:=MinimiseL1ToLinearProgram(abs_coef, constants);
