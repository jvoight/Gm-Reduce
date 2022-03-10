====================================
5T1-[5,5,1]-5-5-11111-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 -1/(x^5 - 1) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    2,
    t - x^5
]
[
    2,
    1,
    -t^2 + x^5
]
[
    3,
    1,
    t^3 - x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 26 

The total time this took was 25.840 


WriteTestToFile(
    filename: 5T1-[5,5,5]-5-5-5-g2.m,
    dextra: 1,
    small_functions_size: 10
)
SortSmallFunctions(
    phi: 1/2/x^3*y + 1/2,
    xs: []
)
ReduceRationalFunction(
    X: Hyperelliptic Curve defined by y^2 = x^6 - 2*x over K,
    phi: 1/2/x^3*y + 1/2,
    P: p
)
Reduction(
    X: Hyperelliptic Curve defined by y^2 = x^6 - 2*x over K,
    p: p
)
Reduction(
    I: Ideal of Graded Polynomial ring of rank 3 over K,
    p: p
)
In file "/opt/magma/current/package/Geometry/Sch/reduction.m", line 39, column 
27:
>>     P := HilbertPolynomial(I);
                             ^
Runtime error in 'HilbertPolynomial': Variable weights are not all 1

====================================
5T2-[5,2,2]-5-221-221-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 1/8/(x^5 - 5/4*x^3 + 5/16*x + 1/16) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    2,
    -25*t*x^4 + 50*t*x^3 - 35*t*x^2 + 10*t*x - t + x^5
]
[
    2,
    7,
    t^2*x^5 - 5*t*x^3 + 5*x - 2
]
[
    3,
    1,
    -5*t^3 + 5*t^2*x^2 - t*x^4 + x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 26 

The total time this took was 17.820 

====================================
5T3-[4,4,2]-41-41-221-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1^2 + 1 over the Rational 
Field defined by
0 

Original map: 
 (1/4107*(-136*nu1 + 4623)*x^5 + 1/37*(-136*nu1 + 4623)*x^4)/(x^5 + (40*nu1 + 
55)*x^4 + (13236*nu1 - 12048)*x^3 + (1006992*nu1 - 709956)*x^2 + (-67346586*nu1 
- 36186777)*x - 7475471046*nu1 - 4016732247) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    1,
    (2*nu1 - 1)*t*x + t + x^5 + (2*nu1 - 1)*x^4
]
[
    2,
    5,
    (2*nu1 - 11)*t^2*x - t^2 + (6*nu1 - 8)*t*x^3 - 2*t*x^2 - x^5 - x^4
]
[
    3,
    1,
    -t^3 + (6*nu1 - 8)*t^2*x^2 - t^2*x + (2*nu1 - 11)*t*x^4 - 2*t*x^3 - x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 30 

The total time this took was 34.030 

====================================
5T3-[5,4,4]-5-41-41-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 - 1/75*x - 118/16875 over Number Field with
defining polynomial $.1^2 + 1 over the Rational Field 

Original map: 
 (256/3125*nu1*x + 1792/46875*nu1)/(x^5 - 1/15*x^4 + 58/225*x^3 + 
1/84375*(-6912*nu1 - 866)*x^2 + 1/6328125*(13824*nu1 - 61343)*x + 
1/2373046875*(960768*nu1 + 309599))*y + (-256/3125*nu1*x^2 + 512/234375*nu1*x + 
1/87890625*(35584*nu1 + 294912))/(x^5 - 1/15*x^4 + 58/225*x^3 + 
1/84375*(-6912*nu1 - 866)*x^2 + 1/6328125*(13824*nu1 - 61343)*x + 
1/2373046875*(960768*nu1 + 309599)) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    2,
    (nu1 + 2)*t^2*x + 2*nu1*t^2 + 8*t*x^5 + 5*nu1*t*x - nu1*t + (11*nu1 - 2)*x^5
]
[
    3,
    1,
    t^3*x^5 - 5*nu1*t^3*x^4 + 40*nu1*t^2*x + (-32*nu1 + 24)*t^2 + 5*t*x + 
        (-3*nu1 - 12)*t + nu1
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 20 

The total time this took was 57.670 

====================================
5T4-[5,2,3]-5-221-311-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 -2/9/(x^5 - 10/3*x^4 + 95/18*x^3 - 25/6*x^2 + 25/16*x - 2/9) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    1,
    -40*t*x^2 + 5*t*x - t - x^5
]
[
    2,
    6,
    3*t^2*x^5 - 10*t*x^3 + 15*x - 8
]
[
    3,
    2,
    -5*t^3*x^5 - 10*t^2*x^3 + t^2*x^2 + 11*t*x + 64
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 26 

The total time this took was 19.060 

====================================
5T4-[5,3,3]-5-311-311-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 2*x^5/(x^5 - 45/8*x^4 + 135/4*x^2 - 729/8) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    4,
    -t*x^5 + 10*x^2 - 15*x + 6
]
[
    2,
    4,
    -t^2*x^5 + t*x^5 - 10*x^2 - 15*x - 36
]
[
    3,
    1,
    15*t^3*x^5 + 30*t^2*x^3 - t^2*x^2 + 7*t*x - 16
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 26 

The total time this took was 19.270 

====================================
5T4-[5,5,2]-5-5-221-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 + 211/375*x - 6214/84375 over Number Field 
with defining polynomial $.1 - 1 over the Rational Field 

Original map: 
 (2048/3125*x - 2048/46875)/(x^5 - 17/15*x^4 + 578/1125*x^3 - 120418/84375*x^2 -
137663/6328125*x - 969431633/2373046875)*y + (-2048/3125*x^2 - 4096/234375*x - 
36800512/87890625)/(x^5 - 17/15*x^4 + 578/1125*x^3 - 120418/84375*x^2 - 
137663/6328125*x - 969431633/2373046875) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    1,
    4*t^2*x^5 - 5*t*x - 2*t + 5*x + 2
]
[
    3,
    1,
    -4*t^2*x^5 + 5*t*x - 2*t - 5*x + 2
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 20 

The total time this took was 14.560 

====================================
5T4-[5,5,3]-5-5-311-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 - 269/375*x + 18962/16875 over Number Field
with defining polynomial $.1 - 1 over the Rational Field 

Original map: 
 (13824/3125*x + 142848/15625)/(x^5 + 47/15*x^4 + 4418/1125*x^3 - 21554/3375*x^2
+ 24529/10125*x - 8693/30375)*y + (-13824/3125*x^2 + 64512/78125*x - 
3896832/390625)/(x^5 + 47/15*x^4 + 4418/1125*x^3 - 21554/3375*x^2 + 
24529/10125*x - 8693/30375) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    7,
    t^2*x^5 - 10*t*x^2 + 15*t*x + 18*t - 24
]
[
    3,
    8,
    3*t^3*x^5 - 10*t^2*x^3 + 15*t*x + 4*t + 4
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 20 

The total time this took was 22.960 


WriteTestToFile(
    filename: 5T4-[5,5,5]-5-5-5-g2.m,
    dextra: 1,
    small_functions_size: 10
)
SortSmallFunctions(
    phi: (x^2 + 243/1024)/(x^5 - 45/32*x^4 + 405/512*x^3 - 3645/16384...,
    xs: []
)
ReduceRationalFunction(
    X: Hyperelliptic Curve defined by y^2 = x^6 + 567/512*x^4 + 111...,
    phi: (x^2 + 243/1024)/(x^5 - 45/32*x^4 + 405/512*x^3 - 3645/16384...,
    P: p
)
Reduction(
    X: Hyperelliptic Curve defined by y^2 = x^6 + 567/512*x^4 + 111...,
    p: p
)
Reduction(
    I: Ideal of Graded Polynomial ring of rank 3 over K,
    p: p
)
In file "/opt/magma/current/package/Geometry/Sch/reduction.m", line 39, column 
27:
>>     P := HilbertPolynomial(I);
                             ^
Runtime error in 'HilbertPolynomial': Variable weights are not all 1

====================================
5T5-[4,2,6]-41-221-32-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 (140625/140608*x^5 - 28125/35152*x^4)/(x^5 - 10/13*x^4 - 335/169*x^3 + 
3420/2197*x^2 + 2160/2197*x - 1728/2197) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    2,
    t*x^2 - 2*t*x + t + 2*x^5 + 25*x^4
]
[
    2,
    4,
    t^2*x^5 - 2*t*x^3 - 25*t*x^2 + x + 27
]
[
    3,
    3,
    -27*t^3 + 25*t^2*x^2 + t^2*x - 2*t*x^3 + x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 30 

The total time this took was 29.480 

====================================
5T5-[4,4,3]-41-41-311-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 (140625/142417*x^5 - 787500/142417*x^4)/(x^5 - 18340/3531*x^4 + 
125440/38841*x^3 + 1756160/142417*x^2 + 3073280/427251*x - 17210368/427251) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    4,
    5*t*x - 3*t + 3*x^5 - 5*x^4
]
[
    2,
    4,
    -8*t^2 - 24000*t*x^4 + 2560*t*x^3 + 64*t + 331776*x^5 + 69120*x^4
]
[
    3,
    3,
    t^3*x^5 - 45*t^2*x^3 + 25*t^2*x^2 + 340*t*x^2 - 210*t*x - 750*x + 486
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 30 

The total time this took was 40.850 

====================================
5T5-[4,6,3]-41-32-311-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1^2 - 6 over the Rational 
Field defined by
0 

Original map: 
 (1/313*(2784*nu1 + 7169)*x^5 + 1/626*(-153*nu1 - 1448)*x^4)/(x^5 + 
1/626*(287*nu1 - 1048)*x^4 + 1/313*(128*nu1 - 282)*x^3 + 1/2817*(-3343*nu1 + 
8172)*x^2 + 1/2817*(-5956*nu1 + 14589)*x + 1/5634*(18283*nu1 - 44784)) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    1,
    (-12*nu1 + 28)*t*x^2 + (-3*nu1 + 7)*t*x + t + x^5 + (nu1 + 1)*x^4
]
[
    2,
    6,
    (-12*nu1 + 28)*t^2*x^2 + (6*nu1 - 19)*t^2*x + 16*t^2 + (-33*nu1 + 81)*t*x^5 
        + 20*t*x^2 + (-9*nu1 - 19)*t*x - x^5
]
[
    3,
    3,
    t^3*x^5 + (-10*nu1 - 30)*t^2*x^3 + (-384*nu1 - 944)*t^2*x^2 + (-4*nu1 + 
        16)*t*x^2 + (22*nu1 + 67)*t*x + (7*nu1 - 13)*x + 32*nu1 + 64
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 30 

The total time this took was 59.810 

====================================
5T5-[5,2,6]-5-2111-32-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 3/2*x^5/(x^5 - 5/2*x^4 - 5/4*x^3 + 45/8*x^2 - 27/8) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    1,
    -t*x^2 + 2*t*x - t + x^5
]
[
    2,
    1,
    -t^2 + 5*t*x^3 - x^5
]
[
    3,
    2,
    -t^3 + 5*t^2*x^2 - x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 30 

The total time this took was 14.300 

====================================
5T5-[5,4,2]-5-41-2111-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 64/63*x^5/(x^5 - 85/504*x^4 - 1445/2268*x^3 - 24565/27216*x^2 + 
1419857/1959552) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    2,
    -t*x^5 + 5*x - 4
]
[
    2,
    4,
    t^2*x^5 - 10*t*x^3 + 25*x - 16
]
[
    3,
    5,
    -t^3*x^5 - 60*t*x^2 + 125*x - 64
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 1 

The number of models that were reduced is 30 

The total time this took was 27.430 

====================================
5T5-[5,4,4]-5-41-41-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 + 1953125/1048576*x + 1220703125/536870912 
over Number Field with defining polynomial $.1 - 1 over the Rational Field 

Original map: 
 (-9765625/1048576*x + 30517578125/1073741824)/(x^5 - 3125/512*x^4 + 
5859375/524288*x^3 - 1220703125/134217728*x^2 + 3814697265625/1099511627776*x - 
286102294921875/562949953421312)*y - 95367431640625/2199023255552/(x^5 - 
3125/512*x^4 + 5859375/524288*x^3 - 1220703125/134217728*x^2 + 
3814697265625/1099511627776*x - 286102294921875/562949953421312) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    1,
    5*t^2*x - 4*t^2 - 2*t*x^5 + x^5
]
[
    3,
    7,
    t^3 - 10*t^2*x^2 - 2*t^2 + 25*t*x^4 + 10*t*x^2 + t + 32*x^5
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 20 

The total time this took was 31.530 

====================================
5T5-[5,6,4]-5-32-41-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 + 1/146484375*(5655936*nu1 - 13761101)*x + 
1/20599365234375*(32784398208*nu1 - 79530479078) over Number Field with defining
polynomial $.1^2 - 6 over the Rational Field 

Original map: 
 (1/30517578125*(-1069645824*nu1 + 2709725184)*x + 
1/95367431640625*(-92066512896*nu1 + 290402574336))/(x^5 + 1/9375*(768*nu1 - 
6913)*x^4 + 1/439453125*(25294848*nu1 - 35320318)*x^3 + 
1/20599365234375*(63299962368*nu1 - 213845723138)*x^2 + 
1/193119049072265625*(27296160448512*nu1 - 54974447094067)*x + 
1/226311385631561279296875*(358703942993689344*nu1 - 992752349979472129))*y + 
(1/30517578125*(1069645824*nu1 - 2709725184)*x^2 + 
1/476837158203125*(288231063552*nu1 - 701587832832)*x + 
1/931322574615478515625*(-31029244749840384*nu1 + 71661730221318144))/(x^5 + 
1/9375*(768*nu1 - 6913)*x^4 + 1/439453125*(25294848*nu1 - 35320318)*x^3 + 
1/20599365234375*(63299962368*nu1 - 213845723138)*x^2 + 
1/193119049072265625*(27296160448512*nu1 - 54974447094067)*x + 
1/226311385631561279296875*(358703942993689344*nu1 - 992752349979472129)) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    1,
    t^2*x^5 + (4*nu1 + 4)*t^2*x^4 + (10*nu1 + 35)*t*x^2 + (-28*nu1 - 73)*t*x + 
        (44*nu1 + 108)*t + 6*nu1 + 9
]
[
    3,
    3,
    t^3*x^5 - 10*t^2*x^3 + (-20*nu1 + 40)*t^2*x^2 + (-90*nu1 - 315)*t^2*x + 
        (40*nu1 + 125)*t*x + (96*nu1 + 240)*t - 32*nu1 - 80
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 20 

The total time this took was 74.540 

====================================
5T5-[5,6,6]-5-32-32-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 - 1953125/459165024*x + 
45166015625/289207845356544 over Number Field with defining polynomial $.1 - 1 
over the Rational Field 

Original map: 
 (-9765625/11019960576*x - 30517578125/1156831381426176)/(x^5 + 3125/52488*x^4 -
37109375/11019960576*x^3 - 28076171875/289207845356544*x^2 + 
41961669921875/7589970693537140736*x - 11539459228515625/1991911908811887214755\
84)*y - 95367431640625/242879062193188503552/(x^5 + 3125/52488*x^4 - 
37109375/11019960576*x^3 - 28076171875/289207845356544*x^2 + 
41961669921875/7589970693537140736*x - 11539459228515625/1991911908811887214755\
84) 

Genus(X) = 1
[ Degree, Index of best model out of 10 (or less), best model]
[
    2,
    1,
    -9*t^2*x^5 + 30*t^2*x^4 - 25*t^2*x^3 + 4*t - 1
]
[
    3,
    3,
    4*t^3*x^5 - 5*t^2*x^2 + 5*t*x^2 - 3*t + 3
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 3 

The number of models that were reduced is 20 

The total time this took was 36.220 

====================================
5T5-[6,6,2]-32-32-221-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 (156250/165997*x^5 - 62500/12769*x^4 + 81250/12769*x^3)/(x^5 - 485/113*x^4 + 
98215/12769*x^3 - 36335/12769*x^2 - 87880/12769*x + 114244/12769) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    9,
    25*t*x^2 + 10*t*x + t - 1024*x^5 - 640*x^4 - 100*x^3
]
[
    2,
    8,
    t^2 - 10*t*x^3 + 5*t*x^2 - 32*x^5
]
[
    3,
    2,
    -12*t^3*x^5 + 15*t^2*x^3 - 5*t*x^2 - 5*t*x + 3
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 30 

The total time this took was 31.690 

====================================
5T5-[6,6,3]-32-32-311-g0.m
Original curve: 
 Curve over Number Field with defining polynomial $.1 - 1 over the Rational 
Field defined by
0 

Original map: 
 (78125/41261*x^5 - 31250/41261*x^4 + 3125/41261*x^3)/(x^5 - 290/341*x^4 + 
815/3751*x^3 - 45/41261*x^2 - 270/41261*x + 27/41261) 

Genus(X) = 0
[ Degree, Index of best model out of 10 (or less), best model]
[
    1,
    2,
    -125*t*x^2 - 25*t*x + t + 16*x^5 + 8*x^4 + x^3
]
[
    2,
    8,
    t^2 - 5*t*x^3 - 5*t*x^2 + x^5
]
[
    3,
    2,
    20*t^3*x^5 - 15*t^2*x^3 - t^2*x^2 + 100*t*x^2 + 3*t*x + 9
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 30 

The total time this took was 38.770 

