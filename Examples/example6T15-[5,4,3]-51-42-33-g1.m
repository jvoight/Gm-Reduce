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

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

SetProfile(true);

ffs:=[];
multiplicities:=[];
support:=[];
for xx in SmallFunctions(PsQsRs, 2*Genus(X)+1) do
  S3orbit:=[ phi, 1/phi, phi-1, 1/(phi-1), 1/phi -1 ];
  for belyimap in [S3orbit[1]] do
    f := ReduceModel(belyimap, xx); f;
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
