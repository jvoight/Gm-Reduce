// Define the base field
K<nu> := RationalsAsNumberField();
// Define the curve
S<x> := PolynomialRing(K);
X := EllipticCurve(S!x^3 - 69295/12582912*x - 15253237/115964116992,S!0);
// Define the map
KX<x,y> := FunctionField(X);
phi := (-600362847/17592186044416*x^2 + 18211006359/36028797018963968*x + 1380901255083/36893488147419103232)/(x^8 - 389/768*x^7 + 321251/4718592*x^6 + 1070293/7247757312*x^5 - 240725312233/712483534798848*x^4 + 16723073889109/1094374709451030528*x^3 + 1527719307251039/6723838214867131564032*x^2 - 152182794452919829/5163907749017957041176576*x + 1341755248500328186849/2030531149437844995903288508416)*y + (600362847/17592186044416*x^4 + 75845839671/9007199254740992*x^3 - 7058065750281/18446744073709551616*x^2 - 145100918243295/4722366482869645213696*x + 236771641668215745/309485009821345068724781056)/(x^8 - 389/768*x^7 + 321251/4718592*x^6 + 1070293/7247757312*x^5 - 240725312233/712483534798848*x^4 + 16723073889109/1094374709451030528*x^3 + 1527719307251039/6723838214867131564032*x^2 - 152182794452919829/5163907749017957041176576*x + 1341755248500328186849/2030531149437844995903288508416);

d := Degree(phi);
phi := 1/phi;
  // should loop over the whole S_3 orbit (1-1/phi, etc.) to find smallest


RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));

PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

//xs := SmallFunctions(X, PsQsRs, 2);
xs := SmallFunctions(PsQsRs, 2);

SetProfile(true);

freds := [];
//for x in SmallFunctions(X, PsQsRs, 2) do
for x in SmallFunctions(PsQsRs, 2) do
    // this is too brutal, surely we want x to have zeros/poles at oo as well;
    // is there any pattern in what works best here?
  f := model(phi, x);
  fred := reducemodel(f);
    // sorry that these are in terms of variables u,v; this helped me to get
    // unconfused about something (fixed a stupid bug),
    // but ultimately it's more meaningful to write this in terms of x,t
  print "f is"; print f; print "f reduced is";  print fred; print "x is"; print x; print "==";
  Append(~freds, <#Sprint(fred), fred, x>);
end for;

ProfilePrintByTotalTime(:Max:=40);

Sort(~freds);
print "Best one:";
print freds[1];

// Possible idea to study: use two small functions (that generate the function field)
// and write phi in terms of these functions, instead of making it a generator


Sort([ <F[1], #Sprint(F[3])> : F in freds ]);
