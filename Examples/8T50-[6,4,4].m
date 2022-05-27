load "config.m";
AttachSpec("../Gm-Reduce/spec");
SetDebugOnError(true);

s := BelyiDBRead("8T50-[6,4,4]-611-422-4211-g0.m");
phi:=BelyiMaps(s)[1];
BestModel(phi);




f:=UnreducedModel(phi);
f:=reducemodel_padic(f);
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
  return [ Log(Abs(Evaluate(x,v : Precision:=prec))) : v in inf_places ];
end function;

mexps := [ Exponents(m) : m in Monomials(f) ];
coefs:=Coefficients(f);
//assert &+[ coefs[i]*(u^mexps[i,1])*v^mexps[i,2] : i in [1..#mexps] ] eq fuv;

UK,mUK:=UnitGroup(K);
k := RealField(prec);
//UU:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) ];
UU:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) and k!0 notin phi(K!(mUK(eps)))  ];
assert UU ne [];

constants := [];
abs_coef := [];
for n in [1..#mexps] do

  alpha_norm := Log(Abs(Norm(coefs[n])))/(r+s);
  log_coef:= phi(coefs[n]);

  for m in [1..r+s] do

    if m le r then
      const:= Log(Abs(Evaluate(coefs[n], inf_places[m]))) - Log(Abs(Norm(coefs[n])))/(r+s);
    else
      const:= Log(Abs(Evaluate(coefs[n], inf_places[m]))) - Log(Abs(Norm(coefs[n])))/(2*(r+s));
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
soln,state:=Solution(L);  //seg faults here


K:=BaseRing(BaseRing(Parent(phi)));
UK, mUK := UnitGroup(Integers(K));
us:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) ];

reducemodel_units_naive := function(f);
  letmeout := false;
  w0 := 1;
  repeat
    S := [];
    oldlen := #Sprint(f);
    for w in us do
      Append(~S, <#Sprint(f*w), w>);
    end for;
    Sort(~S);
    if S[1][1] lt oldlen then
      // print oldlen, S[1];
      oldlen := S[1][1];
      f *:= S[1][2];
      w0 *:= S[1][2];
    else
      letmeout := true;
    end if;
  until letmeout;
  return f, w0;
end function;


intrinsic reducemodel_units_naive(f::RngMPolElt) -> RngMPolElt, SeqEnum
  {Try substituting increasing powers of fundamental units until it stops improving}

  variables:=[ Parent(f).i : i in [1..#Names(Parent(f))] ];

  K:=BaseRing(BaseRing(Parent(phi)));
  UK, mUK := UnitGroup(Integers(K));

  UU := [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) ];
  exp:=1;
  us:= Setseq(Set(&cat[ [ u^e : e in [-exp..exp] ] : u in UU ]));

  yepdone := false;
  oldlen := #Sprint(f);
  oldlen; Sprint(f);
  repeat
    S := [];
    for u,v,w in us do
      tuple:=[u,v,w];
      Append(~S, <#Sprint(Evaluate(f,[tuple[i]*variables[i] : i in [1..#variables]])*tuple[3]), tuple>);
    end for;
    Sort(~S);
    if S[1][1] lt oldlen then
      print oldlen, S[1];
      oldlen := S[1][1];
      exp;
      exp := exp+1;
      old_us:=us;
      us:=Setseq(Set(&cat[ [ u^e : e in [-exp..exp] ] : u in us ]));
      us:= Setseq(Set(us) diff Set(old_us));
    else
      yepdone := true;
    end if;
  until yepdone;

  tuple:=S[1,2];
  return Evaluate(f,[tuple[i]*variables[i] : i in [1..#variables]])*tuple[3]), tuple;
end intrinsic;
