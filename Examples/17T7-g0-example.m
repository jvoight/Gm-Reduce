AttachSpec("spec");
SetDebugOnError(true);
QQ := Rationals();
R<x,y> := PolynomialRing(QQ,2);
F := x*y^17 + 15/4*x*y^16 + 69/10*x*y^15 + 359/40*x*y^14 + 2943/320*x*y^13 + 48301/6400*x*y^12 + 328079/64000*x*y^11 + 745869/256000*x*y^10 + 27691787/20480000*x*y^9 + 41887701/81920000*x*y^8 + 60286491/409600000*x*y^7 + 41159017/1638400000*x*y^6 - 2535961/13107200000*x*y^5 - 602647659/262144000000*x*y^4 - 2137207039/2621440000000*x*y^3 - 1355028301/10485760000000*x*y^2 - 16350628239/1677721600000000*x*y - 1908029761/6710886400000000*x + 8192/263671875*y^2 + 512/48828125;

// integralize F
cs, mons := CoefficientsAndMonomials(F);
dens := [Denominator(el) : el in cs];
printf "F = %o\n\n", F;
F_old := F;
F := LCM(dens)*F;
printf "integralized F = %o\n", F;

F_padic := reducemodel_padic(F);
F_unit := reducemodel_units(F_padic);
