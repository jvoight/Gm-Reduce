
// Belyi maps downloaded from the LMFDB on 21 January 2022.
// Magma code for Belyi map with label 9T32-9_3.3.1.1.1_7.1.1-a

// Geometric data

// Define the base field
R<T> := PolynomialRing(Rationals());
K<nu> := NumberField(R![28, 6, 9, -2, -6, 0, 1]);

// Define the curve
S<x> := PolynomialRing(K);
X := EllipticCurve(S!x^3 + (3354281344/879832127349*nu^5 + 140663240/125690303907*nu^4 - 8537331868/293277375783*nu^3 - 13306304272/879832127349*nu^2 + 3565190068/46306954071*nu + 13594493915/293277375783)*x - 143921036088304/540270595952054289*nu^5 + 8832664077080/77181513707436327*nu^4 + 352794025770268/180090198650684763*nu^3 - 214999362110864/540270595952054289*nu^2 - 147850316828516/28435294523792331*nu + 129693562476298/180090198650684763,S!0);
// Define the map
KX<x,y> := FunctionField(X);
phi := (1/2101140981999124588798317*(1215752832842543139840*nu^5-6617415145255229952*nu^4-9043848325007078792192*nu^3-1907326928863371512832*nu^2+23810188038154585958400*nu+6422149598262798681600)*x^3+1/430076244182454814707361111779*(122549540909780485862328320*nu^5-2710084191444192386865664*nu^4-911155796340064000910709760*nu^3-177284758757214177640506368*nu^2+2399932133316502600895374336*nu+608762594819586736935259648)*x^2+1/264093048578922385974016871660124519*(-17064026718152103635864635661312*nu^5+5910256469367785631111864149504*nu^4+125271745551286525945726070731776*nu^3-16194517240998836472047679683584*nu^2-331589578878095259011442442010624*nu+23484258695630875500925429716480)*x+1/162169241503421659253590774228487722261659*(433924623162942611157328727820935168*nu^5-362686081335711547470531706106240512*nu^4-3127284667574447317566242563127553024*nu^3+1974769246585100754227913443127387136*nu^2+8342334473227841462760776974345753600*nu-4728233829574181466808979920247600640))/(x^9+1/204687*(688*nu^5+2920*nu^4+4*nu^3-12592*nu^2-4636*nu-26101)*x^8+1/15435651357*(142243536*nu^5-43391152*nu^4-1108418240*nu^3-52163904*nu^2+2929665328*nu+692625380)*x^7+1/25727171235812109*(-33967112877904*nu^5+114663965536400*nu^4+233752310058096*nu^3-760274006642816*nu^2-660025104923456*nu+1906985732270244)*x^6+1/5266017498744673154883*(-2251212913491365776*nu^5-2374562502856766992*nu^4+17005397510997585920*nu^3+20206920876971358784*nu^2-44026961439098157872*nu-55123242498963680726)*x^5+1/1077885323765550914053536621*(50292428283252130772368*nu^5-32198629388515484809600*nu^4-362637331762051006889496*nu^3+165680399508567340522976*nu^2+978585570914507727693112*nu-360395007020893023730626)*x^4+1/661887339796797959834628751027881*(44764623807884296723721136*nu^5+2145435772890242133292295472*nu^4-764654693537153584406017472*nu^3-15706158317553793406033833152*nu^2+668718139229498164491159504*nu+40044253130770654747066391204)*x^3+1/45159911306995728001556885053881292749*(-3511028203345554045622771143280*nu^5-1457654366106606707285626603664*nu^4+25870276384177203760363774417072*nu^3+15213808512883214342461758382080*nu^2-66778060838325957840424692270176*nu-41462862495254472176732371133612)*x^2+1/27730940297085103732364022393071400506743689*(28601135721701912788833783681455984*nu^5-34915214335735421581684199869792528*nu^4-198283316885110982856812148237676672*nu^3+222766485955999581227793881399133760*nu^2+537642537946592458994914902784714576*nu-556000666141015540682268831880917887)*x+1/119199422508378631180994287682982720866000754877203*(457199704474074388313562029968887085024*nu^5+1346519342883608867157075304166163804856*nu^4-3727295422130447876233598545058348624876*nu^3-10689541235678022556864415988622033474096*nu^2+8869363651367450195678772716266834185284*nu+27743997280443495905237200957039674829207))*y+(1/2101140981999124588798317*(-1215752832842543139840*nu^5+6617415145255229952*nu^4+9043848325007078792192*nu^3+1907326928863371512832*nu^2-23810188038154585958400*nu-6422149598262798681600)*x^4+1/10506148250742824759279821444887*(-209838183766482101586829312*nu^5-966181488053784551045797888*nu^4+1831407285208740442454511616*nu^3+7455932489418661872873017344*nu^2-4540461401170053886756032512*nu-19969906227941122684269082624)*x^3+1/5017767922999525333506320561542365861*(233443163056862518819092658337792*nu^5+282922030879002332830788102720512*nu^4-1814555446660887807891975742486528*nu^3-2456929120276229493303824872069120*nu^2+4694649173843258929726195072004096*nu+6763697394555808182066262738439168)*x^2+1/3081215588565011525818224710341266722971521*(-24437877452726805586420107202602835968*nu^5-12232095366599619321430909776814688256*nu^4+185155150318725431302905745935249010688*nu^3+129218939660581412015369653170897825792*nu^2-483965366480860486016131621783849660416*nu-369608640187620740121731346540325480448)*x+1/13244380278708736797888254186998080096222306097467*(4170054086255216918027777305446819473591296*nu^5-214373223561901923033023532072364258474496*nu^4-30958570307772740976563898994990260472699904*nu^3-5105022167189556169832646688547214889028608*nu^2+81589919856216023263246789943464844478654464*nu+18258811878921224538565358928411638363366912))/(x^9+1/204687*(688*nu^5+2920*nu^4+4*nu^3-12592*nu^2-4636*nu-26101)*x^8+1/15435651357*(142243536*nu^5-43391152*nu^4-1108418240*nu^3-52163904*nu^2+2929665328*nu+692625380)*x^7+1/25727171235812109*(-33967112877904*nu^5+114663965536400*nu^4+233752310058096*nu^3-760274006642816*nu^2-660025104923456*nu+1906985732270244)*x^6+1/5266017498744673154883*(-2251212913491365776*nu^5-2374562502856766992*nu^4+17005397510997585920*nu^3+20206920876971358784*nu^2-44026961439098157872*nu-55123242498963680726)*x^5+1/1077885323765550914053536621*(50292428283252130772368*nu^5-32198629388515484809600*nu^4-362637331762051006889496*nu^3+165680399508567340522976*nu^2+978585570914507727693112*nu-360395007020893023730626)*x^4+1/661887339796797959834628751027881*(44764623807884296723721136*nu^5+2145435772890242133292295472*nu^4-764654693537153584406017472*nu^3-15706158317553793406033833152*nu^2+668718139229498164491159504*nu+40044253130770654747066391204)*x^3+1/45159911306995728001556885053881292749*(-3511028203345554045622771143280*nu^5-1457654366106606707285626603664*nu^4+25870276384177203760363774417072*nu^3+15213808512883214342461758382080*nu^2-66778060838325957840424692270176*nu-41462862495254472176732371133612)*x^2+1/27730940297085103732364022393071400506743689*(28601135721701912788833783681455984*nu^5-34915214335735421581684199869792528*nu^4-198283316885110982856812148237676672*nu^3+222766485955999581227793881399133760*nu^2+537642537946592458994914902784714576*nu-556000666141015540682268831880917887)*x+1/119199422508378631180994287682982720866000754877203*(457199704474074388313562029968887085024*nu^5+1346519342883608867157075304166163804856*nu^4-3727295422130447876233598545058348624876*nu^3-10689541235678022556864415988622033474096*nu^2+8869363651367450195678772716266834185284*nu+27743997280443495905237200957039674829207));

print "phi has divisor";
Support(Divisor(phi));

load "../../reducecurve.m";

RsandPs := Support(Divisor(phi));
RsandQs := Support(Divisor(phi-1));
PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
printf "ramification points = %o\n", PsQsRs;

print "computing small functions supported at points above";
xs := SmallFunctions(PsQsRs, 2);
for x_op in xs do
  pts, mults := Support(Divisor(x_op));
  printf "x_op has support\n%o,\n %o\n", pts, mults;
  time F_res := PlaneModel(phi, x_op);
  print F_res;
  print "-------------------------------------";
end for;