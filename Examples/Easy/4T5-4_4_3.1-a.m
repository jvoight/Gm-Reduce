AttachSpec("spec");
// Belyi maps downloaded from the LMFDB on 28 January 2022.
// Magma code for Belyi map with label 4T5-4_4_3.1-a


// Group theoretic data

d := 4;
i := 5;
G := TransitiveGroup(d,i);
sigmas := [[Sym(4) | [2, 3, 4, 1], [3, 1, 4, 2], [3, 2, 4, 1]]];
embeddings := [ComplexField(15)![1.0, 0.0]];

// Geometric data

// Define the base field
K<nu> := RationalsAsNumberField();
// Define the curve
S<x> := PolynomialRing(K);
X := EllipticCurve(S!x^3 + 47/768*x + 2359/55296,S!0);
// Define the map
KX<x,y> := FunctionField(X);
phi := 27/64/(x^4-13/12*x^3-155/384*x^2-1225/27648*x-8375/5308416)*y+(-27/64*x^2+9/512*x+1401/16384)/(x^4-13/12*x^3-155/384*x^2-1225/27648*x-8375/5308416);



RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

SetProfile(true);
small_functions:=SmallFunctions(PsQsRs, 2*Genus(X)+1);

ffs:=[];
multiplicities:=[];
support:=[];
for xx in SmallFunctions(PsQsRs, 2*Genus(X)+1) do
  S3orbit:=[ phi, 1/phi, 1-phi, phi/(phi-1), 1-1/phi, 1/(1-phi)  ];
  for belyimap in [S3orbit[1]] do
    f := model(belyimap, xx); f;
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
