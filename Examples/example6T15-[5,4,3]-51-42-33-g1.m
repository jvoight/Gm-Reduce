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


fuv:=1/95551488*(-11132812500*nu + 35205078125)*t^2 + 1/3456*(884375*nu -
    2796875)*t*x^3 + 1/110592*(-136250000*nu + 430859375)*t*x^2 +
    1/442368*(975859375*nu - 3085937500)*t*x + 1/23887872*(-30630859375*nu +
    96863281250)*t + x^6 + 1/4*(-15*nu + 50)*x^5;

reducemodel_padic(fuv);
(-1824000000000*nu + 5768000000000)*t^2 + (2445120000*nu - 7732800000)*t*x^3 +
    (3722625000*nu - 11772000000)*t*x^2 + (2107856250*nu - 6665625000)*t*x +
    (387453125*nu - 1225234375)*t + 5832*x^6 + (7290*nu - 21870)*x^5

reducemodel_padic_old(fuv);
    (-3648*nu + 11536)*t^2 + (8824*nu - 27904)*t*x^3 + (-19750*nu + 62455)*t*x^2 +
        (16441*nu - 51991)*t*x + (-4443*nu + 14050)*t + (-6*nu + 19)*x^6 + (-68*nu +
        215)*x^5


        > [ [ Valuation(cc,pp) : cc in Coefficients(reducemodel_padic(fuv)) ] : pp in SS ];
        [
            [ 24, 18, 7, 3, 0, 6, 2 ],
            [ 0, 3, 3, 3, 3, 6, 6 ],
            [ 0, 6, 6, 7, 0, 6, 6 ],
            [ 18, 9, 12, 11, 13, 0, 2 ]
        ]
        > [ [ Valuation(cc,pp) : cc in Coefficients(reducemodel_padic_old(fuv)) ] : pp\
         in SS ];
        [
            [ 8, 7, 0, 0, 1, 0, 0 ],
            [ 0, 0, 0, 0, 0, 0, 0 ],
            [ 0, 3, 4, 6, 0, 0, 1 ],
            [ 0, 0, 2, 0, 1, 0, 1 ]
        ]
       reducemodel_padic_old(fuv : Integral:=false, Polyhedron:=false);
        1/81*(-5330*nu + 16855)*t^2 + (-38*nu + 120)*t*x^3 + 1/4*(-2435*nu + 7700)*t*x^2
    + 1/2*(-7259*nu + 22955)*t*x + (-7025*nu + 22215)*t + 1/9*(4*nu + 13)*x^6 +
    nu*x^5
