AttachSpec("code/spec_database");
load "../../Gm Reduce/Code/reducecurve.m";

filename:="6T15-[5,4,3]-51-42-33-g1.m";
X:=BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
phi:=BelyiDBRead(filename)`BelyiDBBelyiMaps[1];
X; phi;

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

SetProfile(true);

ffs:=[];
multiplicities:=[];
support:=[];
for xx in SmallFunctions(PsQsRs, 2*Genus(X)+1) do
  S3orbit:=[ phi, 1/phi, phi-1, 1/(phi-1), 1/phi -1 ];
  for belyimap in S3orbit do
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
