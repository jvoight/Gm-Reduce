//printf "%m \n", K;
AttachSpec("spec");
SetDebugOnError(true);
K<nu>:=ext<K|Polynomial(K, [-2,0,1])> where K is RationalField();
Rtx<t,x>:= PolynomialRing(K,2);

UK, mUK := UnitGroup(Integers(K));
eps:= K!mUK(UK.2);

f:=5*x^3-2*eps^7*t*x;

fu, tup:=reducemodel_units(f : prec:=10);
//tup is huge

prec:=10;
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

constants := [];
abs_coef := [];
for n in [1..#mexps] do

  alpha_norm := Log(k!Abs(Norm(coefs[n])))/(r+s);
  log_coef:= phi(coefs[n]);

  for m in [1..r+s] do

    if m le r then
      const:= Log(k!Abs(Evaluate(coefs[n], inf_places[m]))) - Log(k!Abs(Norm(coefs[n])))/(r+s);
    else
      const:= Log(k!Abs(Evaluate(coefs[n], inf_places[m]))) - Log(k!Abs(Norm(coefs[n])))/(2*(r+s));
    end if;
    Append(~constants, [const]);

    etas:= [ phi(eps)[m] : eps in UU ];
    lhs:=&cat[ [ eta*a : a in mexps[n] ] cat [eta] : eta in etas ];
    Append(~abs_coef, lhs);

  end for;
end for;

constants:=Matrix(k,constants);
abs_coef:=Matrix(k,abs_coef);
L:=MinimiseL1ToLinearProgram(abs_coef, constants);
soln,state:=Solution(L);

obj:=ObjectiveFunction(L);
EvaluateAt(L,soln);
n:=2;
EvaluateAt(L,Parent(soln)![6*n,3*n,-9*n,0,0,0,0]);
EvaluateAt(L,Parent(soln)![6*n,3*n,-9*n,0,0,0,0]);


rows:=[];
nos:=[];
for i in [1..NumberOfConstraints(L)] do
  row,no:=Constraint(L,i);
  Append(~rows,Eltseq(row));
  Append(~nos,Eltseq(no));
end for;

PolyhedronWithInequalities(rows,nos);
