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

lc:=Coefficients(g)[24];
Factorization(lc*Integers(K)); //[ p2,1^23, p2,2^2, p2,3^23 ]
fac:= [a[1] : a in Factorization(lc*Integers(K)) ];

/*for these primes the minimal points are
    [        (-3, -1, 0)    ],
    [        (-2, 0, -2),   (-14, -1, 21)   ],
    [        (0, -1, 0)    ]]
    and the objective function at these points are [-345,-94,-94,-276]
    Note absorbing the leading coefficient is equivalent to [23,0,-23] and the objective function
    is only -23 at this point. So the LP is finding a much better solution, but we still want to clear denominators!!*/
[ CoefficientValuations(g)[i] : i in [1..3] ];
/*[ 25, 24, 23, 23, 21, 20, 19, 21, 18, 17, 16, 16, 14, 13, 12, 13, 9, 8, 7, 7, 5, 4, 3, 23 ],
[ 19, 14, 14, 12, 13, 14, 15, 12, 14, 11, 11, 9, 10, 11, 11, 8, 14, 6, 6, 4, 5, 6, 7, 2 ],
[ 38, 34, 28, 27, 23, 21, 16, 18, 16, 17, 13, 14, 14, 12, 8, 10, 10, 8, 5, 7, 4, 4, 0, 23 ]*/

/*After finding small elements in the ideals at just these primes we get
new_fuv:=
(132639266537472*nu1^3 + 117930443374592*nu1^2 + 1090872379736064*nu1 -
    2008080488726528)*t*x^22 + (-928266836622336*nu1^3 - 938605281375232*nu1^2 -
    7343720024024064*nu1 + 6828594803814400)*t*x^21 + (1079464369231872*nu1^3 +
    3335850303373056*nu1^2 + 15270998090936832*nu1 - 7312689831585792)*t*x^20 +
    (-3020016195105024*nu1^3 - 1780180196029952*nu1^2 - 17388054250930944*nu1 -
    603852848475136)*t*x^19 + (-1962674489611648*nu1^3 + 4198689384492704*nu1^2
    + 6777297043923520*nu1 + 8504330415460096)*t*x^18 + (-2775554532333408*nu1^3
    + 1271758516305552*nu1^2 + 347480071484448*nu1 - 11891353533393024)*t*x^17 +
    (-3184090483300980*nu1^3 + 99097405359859*nu1^2 - 3448052908810866*nu1 +
    4970847387150728)*t*x^16 + (-1405199326131440*nu1^3 + 66438634006680*nu1^2 +
    4324579058953520*nu1 - 4391285105452224)*t*x^15 + (-1006048963369712*nu1^3 -
    893329399326372*nu1^2 - 241544775246504*nu1 - 725123429357152)*t*x^14 +
    (-296740601646816*nu1^3 - 383912163113136*nu1^2 + 1012211341086240*nu1 -
    543476926164096)*t*x^13 + (-96792553688040*nu1^3 - 226782032477202*nu1^2 +
    235352116179948*nu1 - 332457134593584)*t*x^12 + (-21690183251376*nu1^3 -
    71607260632344*nu1^2 + 78131095316784*nu1 - 37167887870400)*t*x^11 +
    (-2289884611776*nu1^3 - 17175614288720*nu1^2 + 27077390628384*nu1 -
    23572931943040)*t*x^10 + (-415151582768*nu1^3 - 3666765818952*nu1^2 +
    2345002505456*nu1 - 1846892812224)*t*x^9 + (28207380228*nu1^3 -
    449246818119*nu1^2 + 769944621162*nu1 - 398264048232)*t*x^8 +
    (3947239392*nu1^3 - 58182190120*nu1^2 + 37502475936*nu1 - 50903626496)*t*x^7
    + (844212608*nu1^3 - 4092482736*nu1^2 + 4293364000*nu1 + 806070144)*t*x^6 +
    (115387044*nu1^3 - 218047452*nu1^2 + 422528124*nu1 - 461314128)*t*x^5 +
    (3182004*nu1^3 - 12000158*nu1^2 - 18248016*nu1 + 19880096)*t*x^4 +
    (412137*nu1^3 + 157987*nu1^2 + 1476255*nu1 - 394036)*t*x^3 + (3450*nu1^3 -
    5382*nu1^2 - 51474*nu1 - 1656)*t*x^2 + (138*nu1^3 + 322*nu1^2 + 138*nu1 +
    1196)*t*x + 1/2*(3*nu1^3 + 5*nu1^2 + 3*nu1 - 26)*t + (-2382*nu1^3 -
    2227*nu1^2 - 19097*nu1 + 29196)*x^23;
Still has a leading coefficient!*/
[ CoefficientValuations(new_fuv)[i] : i in [1..3] ];
/*[ 44, 42, 40, 39, 36, 34, 32, 33, 29, 27, 25, 24, 21, 19, 17, 17, 12, 10, 8,7, 4, 2, 0, 44 ],
[ 16, 11, 11, 9, 10, 11, 12, 9, 11, 8, 8, 6, 7, 8, 8, 5, 11, 3, 3, 1, 2, 3, 4, 1 ],
[ 16, 13, 8, 8, 5, 4, 0, 3, 2, 4, 1, 3, 4, 3, 0, 3, 4, 3, 1, 4, 2, 3, 0, 0 ]
The first prime has gotten worse!*/

fu1:=reducemodel_units(new_fuv);
fp2:=reducemodel_padic(fu1: Integral:=true,ClearDenominators:=false,Speedy:=false );
fu2:=reducemodel_units(fp2);
fp3:=reducemodel_padic(fu2: Integral:=true,ClearDenominators:=false,Speedy:=false );
fu3:=reducemodel_units(fp3);

adding the condition that x=0 after reducing gives us something with leading coefficient which is a unit, but reducing by units doesn't get rid of it!



clearing leading coef of fu3 and factoring  gives
(1)*(x)^23 + (1/452984832*(-3093833761*nu1^3 - 3181321403*nu1^2 -
25681700327*nu1 + 34940538644))*(1/6*(-19*nu1^3 - 35*nu1^2 - 179*nu1 + 38)*x^2 +
1/3*(2*nu1^3 + 25*nu1^2 - 5*nu1 + 8)*x + 1/6*(5*nu1^3 - 5*nu1^2 + nu1 +
2))^2*(1/6*(-31*nu1^3 - 17*nu1^2 - 305*nu1 + 182)*x^2 + 1/2*(5*nu1^3 + 9*nu1^2 +
13*nu1 + 6)*x + 1/6*(nu1^3 - 7*nu1^2 - nu1 - 2))^1*(1/6*(-77*nu1^3 - 193*nu1^2 -
601*nu1 + 694)*x^4 + 1/3*(-34*nu1^3 + 52*nu1^2 + 133*nu1 - 322)*x^3 +
1/3*(20*nu1^3 + 181*nu1^2 - 209*nu1 + 230)*x^2 + (8*nu1^3 - 12*nu1^2 + 7*nu1 -
6)*x + 1/2*(-nu1^3 - 7*nu1^2 + 9*nu1 - 6))^4*t

(1)*(x)^23 + (1/1296*(74113*nu1^3 + 81117*nu1^2 + 625509*nu1 -
785230))*((-17*nu1^3 - 39*nu1^2 - 187*nu1 - 32)*x^2 + 1/6*(11*nu1^3 + 13*nu1^2 +
145*nu1 - 28)*x + 1/6*(-nu1^3 + nu1^2 - 5*nu1 + 2))^2*(1/3*(26*nu1^3 + 58*nu1^2
+ 238*nu1 - 196)*x^2 + 1/6*(-nu1^3 + nu1^2 - 83*nu1 + 32)*x + 1/6*(nu1^3 - nu1^2
+ 5*nu1 - 2))^1*(1/3*(254*nu1^3 + 682*nu1^2 + 3154*nu1 + 4616)*x^4 + (-108*nu1^3
- 220*nu1^2 - 636*nu1 - 272)*x^3 + 1/3*(14*nu1^3 + 310*nu1^2 + 154*nu1 + 56)*x^2
+ (4*nu1^3 + 2*nu1^2 - 22*nu1 + 16)*x + 1/6*(nu1^3 - 7*nu1^2 + 11*nu1 - 8))^4*t
