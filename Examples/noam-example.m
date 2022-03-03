//AttachSpec("code/spec_database");
//load "../Gm-Reduce/reducecurve.m";
AttachSpec("spec");

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
//xs := SmallFunctions(PsQsRs, 2*Genus(X)+1);
//"The number of small functions is"; #xs;
small_functions:=SmallFunctions(PsQsRs, 1);

SetProfile(true);
ffs:=[];
multiplicities:=[];
support:=[];
for xx in small_functions do
  // this is too brutal, surely we want x to have zeros/poles at oo as well;
  // is there any pattern in what works best here?
  //try 1/phi etc
  S3orbit:=[ phi, 1/phi, 1-phi, phi/(phi-1), 1-1/phi, 1/(1-phi)  ];
  for belyimap in S3orbit do
    f := model(belyimap, xx);

    Append(~ffs,<#Sprint(f),f, Index(S3orbit, belyimap), Degree(xx) >);
    sup,mult:=Support(Divisor(xx));
    Append(~support, mult);
    //insert reduction into model
  end for;
end for;
ParallelSort(~ffs,~support);

shortest_ffs:=[ ffs[i] : i in [1..Min(#ffs,10)] ];
for fuv_es in shortest_ffs do
  fuv := fuv_es[2];

  padic_redfuv:= reducemodel_padic(fuv);
  unit_redfuv := reducemodel_units(padic_redfuv);
  fuv_display := PolynomialToFactoredString(MultivariateToUnivariate(unit_redfuv));
  print fuv_display;
  print "";
end for;

ProfilePrintByTotalTime(:Max:=40);





fuv:=shortest_ffs[1,2];
fp1:=reducemodel_padic(fuv : Integral:=false, ClearDenominators:=false); //fairly quick
fp2:=reducemodel_padic(fuv : Integral:=true, ClearDenominators:=true);  //slow due to polyhedra, in particular meet, reduct, MinimalRGenerators(), Polyhedron(), Inequalities took a long time.
fp3:=reducemodel_padic(fuv : Integral:=true, ClearDenominators:=false); //~2 mins, same as fp2 but even longer because its 3-dimensional

CoefficientValuations(fp1);
CoefficientValuations(fp2);
CoefficientValuations(fp3);

fu1:=reducemodel_units(fp1);
fu2:=reducemodel_units(fp2);
fu3:=reducemodel_units(fp3);

K<nu1>:=F;
_<t,x>:=PolynomialRing(K,2);

//why is it not clearing denominators??
g:=(-352330873874*nu1^3 - 140724214886*nu1^2 - 2459611303550*nu1 +
6316249625192)*(x)^23 + (1/324*(-176407*nu1^3 - 117699*nu1^2 - 1330647*nu1 +
2664154))*((17*nu1^3 + 39*nu1^2 + 187*nu1 + 32)*x^2 + 1/6*(-11*nu1^3 - 13*nu1^2
- 145*nu1 + 28)*x + 1/6*(nu1^3 - nu1^2 + 5*nu1 - 2))^2*(1/3*(26*nu1^3 + 58*nu1^2
+ 238*nu1 - 196)*x^2 + 1/6*(-nu1^3 + nu1^2 - 83*nu1 + 32)*x + 1/6*(nu1^3 - nu1^2
+ 5*nu1 - 2))^1*(1/3*(254*nu1^3 + 682*nu1^2 + 3154*nu1 + 4616)*x^4 + (-108*nu1^3
- 220*nu1^2 - 636*nu1 - 272)*x^3 + 1/3*(14*nu1^3 + 310*nu1^2 + 154*nu1 + 56)*x^2
+ (4*nu1^3 + 2*nu1^2 - 22*nu1 + 16)*x + 1/6*(nu1^3 - 7*nu1^2 + 11*nu1 - 8))^4*t;


fp1:=reducemodel_padic(g);
lc:=Coefficients(fp1)[24];
assert IsUnit(lc);
fu1:=reducemodel_units(fp1);
assert Coefficients(fu1)[24]^2 eq 1;
