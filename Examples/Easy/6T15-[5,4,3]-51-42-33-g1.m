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

BestModel(phi);

/*
RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

SetProfile(true);

ffs:=[];
multiplicities:=[];
support:=[];

effort:=15; degree:=2;
//Ceiling((Genus(x)+3)/2);
//that's a minimum, keep going if it's not good/ there aren't enough functions?
xs:=SmallFunctions(PsQsRs, degree);
xs_sorted:=SortSmallFunctions(phi,xs);
#xs;

reduced_models:=[];
for xx in [ xs_sorted[i] : i in [1..#xs] ] do
  fred:=ReduceModel(phi, xx);
  Append(~reduced_models,<#Sprint(fred),fred>);
end for;

reduce_model_sorted:=Sort(reduced_models);
f:=reduce_model_sorted[1];
Index(reduced_models,f);
reduced_models[1,2];


for xx in SmallFunctions(PsQsRs, 2*Genus(X)+1) do
  S3orbit:=[ phi, 1/phi, phi-1, 1/(phi-1), 1/phi -1 ];
  for belyimap in [S3orbit[1]] do
    f := model(belyimap, xx);
    Append(~ffs,<#Sprint(f),f, Index(S3orbit, belyimap), Degree(xx) >);
    sup,mult:=Support(Divisor(xx));
    Append(~support, mult);
  end for;
end for;

ParallelSort(~ffs,~support);
shortest_ffs:=[ ffs[i] : i in [1..Min(#ffs,5)] ];
fuv:=shortest_ffs[1,2];

padic_redfuv:= reducemodel_padic(fuv);
unit_redfuv := reducemodel_units(padic_redfuv);
fuv_display := PolynomialToFactoredString(MultivariateToUnivariate(unit_redfuv));

ProfilePrintByTotalTime(:Max:=40);

_<t,x>:=PolynomialRing(K,2);
fuv:=1/95551488*(-11132812500*nu + 35205078125)*t^2 + 1/3456*(884375*nu -
    2796875)*t*x^3 + 1/110592*(-136250000*nu + 430859375)*t*x^2 +
    1/442368*(975859375*nu - 3085937500)*t*x + 1/23887872*(-30630859375*nu +
    96863281250)*t + x^6 + 1/4*(-15*nu + 50)*x^5;

fp1:=reducemodel_padic(fuv);
fu1:=reducemodel_units(fp1);
/* fu1:= 16*t^2 + (40*nu - 112)*t*x^3 + (-10*nu + 55)*t*x^2 + (-8*nu - 37)*t*x + (3*nu +
    10)*t + x^6 + (nu + 5)*x^5 */
*/