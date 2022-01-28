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
    f := model(belyimap, xx); f;

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
