
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


WriteTestToFile(
    filename: 6T10-[4,4,2]-42-42-2211-g0.m,
    dextra: 3,
    small_functions_size: 20
)
ReducedModelS3Orbit(
    phi: (x^6 - 2*x^5 + x^4)/(x^6 - 2*x^5 + x^4 + 1/9*(2*nu1 + 3)*x^2...,
    x_op: (x + 1/9*(-2*nu1 - 3))/(x^3 - x^2)
)
ReducedModelsS3Orbit(
    phi: (x^6 - 2*x^5 + x^4)/(x^6 - 2*x^5 + x^4 + 1/9*(2*nu1 + 3)*x^2...,
    x_op: (x + 1/9*(-2*nu1 - 3))/(x^3 - x^2)
)
S3Orbit(
    f: (21*nu1 - 36)*t^3*x^6 + (-9*nu1 + 18)*t^3*x^4 + 3*nu1*t^3*x^...
)
S3Action(
    tau: Id($),
    f: (21*nu1 - 36)*t^3*x^6 + (-9*nu1 + 18)*t^3*x^4 + 3*nu1*t^3*x^...
)
ComputeThirdRamificationValue(
    f: (21*nu1 - 36)*t^3*x^6 + (-9*nu1 + 18)*t^3*x^4 + 3*nu1*t^3*x^...
)
In file "/scratch/home/cs/Gm-Reduce/reducecurve.m", line 140, column 27:
>>   KC<t,x> := FunctionField(C);
                             ^
Runtime error in 'FunctionField': Scheme must be defined by an irreducible 
polynomial


WriteTestToFile(
    filename: 6T10-[4,4,3]-42-42-3111-g0.m,
    dextra: 3,
    small_functions_size: 20
)
ReducedModelS3Orbit(
    phi: (6561/161*x^6 - 74358/161*x^5 + 210681/161*x^4)/(x^6 - 9078/...,
    x_op: (x^3 - 867/161*x^2 - 867/161*x + 4913/161)/(x^3 - 17/3*x^2)
)
ReducedModelsS3Orbit(
    phi: (6561/161*x^6 - 74358/161*x^5 + 210681/161*x^4)/(x^6 - 9078/...,
    x_op: (x^3 - 867/161*x^2 - 867/161*x + 4913/161)/(x^3 - 17/3*x^2)
)
S3Orbit(
    f: -t^3*x^6 - 3*t^3*x^5 - 3*t^3*x^4 - t^3*x^3 - 3*t^2*x^4 - 6*t...
)
S3Action(
    tau: Id($),
    f: -t^3*x^6 - 3*t^3*x^5 - 3*t^3*x^4 - t^3*x^3 - 3*t^2*x^4 - 6*t...
)
ComputeThirdRamificationValue(
    f: -t^3*x^6 - 3*t^3*x^5 - 3*t^3*x^4 - t^3*x^3 - 3*t^2*x^4 - 6*t...
)
In file "/scratch/home/cs/Gm-Reduce/reducecurve.m", line 140, column 27:
>>   KC<t,x> := FunctionField(C);
                             ^
Runtime error in 'FunctionField': Scheme must be defined by an irreducible 
polynomial

====================================
6T10-[4,4,3]-42-42-33-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 + 1007401764096/602425897921*x + 
273179726831591424/467579487356261281 over Number Field with defining polynomial
$.1 - 1 over the Rational Field 

Original map: 
 (1077632058507264/530737216068401*x^2 - 
472054567964989980672/411937528360866188561*x - 
224014439336921925307858944/319729843950098261779694321)/(x^6 - 
1314144/776161*x^5 + 3022205292288/602425897921*x^4 - 
2227465464934514688/467579487356261281*x^3 + 
2464655906166208312836096/362916962485923112122241*x^2 - 
873992335073000891588611670016/281681992520036568627910696801*x + 
542371675391582017291361566838489088/218630576996344103142807794339760961)*y + 
(1726974452736/602425897921*x^4 - 1344884809017065472/467579487356261281*x^3 + 
2595828069239479903715328/362916962485923112122241*x^2 - 
1403226448006478940128428425216/281681992520036568627910696801*x + 
419060223852319425531855995871952896/218630576996344103142807794339760961)/(x^6 
- 1314144/776161*x^5 + 3022205292288/602425897921*x^4 - 
2227465464934514688/467579487356261281*x^3 + 
2464655906166208312836096/362916962485923112122241*x^2 - 
873992335073000891588611670016/281681992520036568627910696801*x + 
542371675391582017291361566838489088/218630576996344103142807794339760961) 

Genus(X) = 1
[ Degree, Index of best model out of 20 (or less), best model]
[
    2,
    1,
    t^2 - 8*t + x^6 - 6*x^5 + 9*x^4
]
[
    3,
    2,
    t + x^2
]
[
    4,
    2,
    t + x^2
]
[
    5,
    2,
    t^5*x^6 - 6*t^4*x^5 + 9*t^3*x^4 + t - 8
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 3 

The number of models that were reduced is 80 

The total time this took was 6666.170 


WriteTestToFile(
    filename: 6T11-[6,2,4]-6-2211-411-g0.m,
    dextra: 3,
    small_functions_size: 20
)
ReducedModelS3Orbit(
    phi: 1/3*x^6/(x^2 + 2/3),
    x_op: x^2
)
ReducedModelsS3Orbit(
    phi: 1/3*x^6/(x^2 + 2/3),
    x_op: x^2
)
S3Orbit(
    f: t^2*x^2 + 2*t^2*x + t^2 - 2*t*x^4 - 2*t*x^3 + x^6
)
S3Action(
    tau: Id($),
    f: t^2*x^2 + 2*t^2*x + t^2 - 2*t*x^4 - 2*t*x^3 + x^6
)
ComputeThirdRamificationValue(
    f: t^2*x^2 + 2*t^2*x + t^2 - 2*t*x^4 - 2*t*x^3 + x^6
)
In file "/scratch/home/cs/Gm-Reduce/reducecurve.m", line 140, column 27:
>>   KC<t,x> := FunctionField(C);
                             ^
Runtime error in 'FunctionField': Scheme must be defined by an irreducible 
polynomial

====================================
6T11-[6,4,2]-6-42-222-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 - 52/27*x - 560/729 over Number Field with 
defining polynomial $.1 - 1 over the Rational Field 

Original map: 
 32/27/(x^3 - 2/3*x^2 - 32/27*x + 640/729) 

Genus(X) = 1
[ Degree, Index of best model out of 20 (or less), best model]
[
    2,
    1,
    t*x + t - x^3
]
[
    3,
    1,
    2*t^3 - 3*t^2*x^2 - 6*t^2 + 9*t*x^2 + x^6
]
[
    4,
    1,
    t^2 - t*x^2 - x^3
]
[
    5,
    17,
    27*t^5*x^6 + 54*t^4*x^4 + 27*t^2*x^2 - 2*t + 2
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 76 

The total time this took was 698.200 

====================================
6T11-[6,4,4]-6-411-42-g1.m
Original curve: 
 Elliptic Curve defined by y^2 = x^3 - 175/432*x - 625/11664 over Number Field 
with defining polynomial $.1 - 1 over the Rational Field 

Original map: 
 125/432/(x^3 - 5/6*x^2 - 125/432*x + 3125/11664) 

Genus(X) = 1
[ Degree, Index of best model out of 20 (or less), best model]
[
    2,
    1,
    t*x + t - x^3
]
[
    3,
    2,
    t^3*x^6 - 3*t^2*x^4 - 45*t^2*x^2 + 3*t*x^2 + 15*t - 1
]
[
    4,
    2,
    t^2 - t*x^2 - x^3
]
[
    5,
    2,
    5*t^5*x^6 - 3*t^4*x^4 + 12*t - 225
]
Ceiling((Genus(X)+3)/2) = 2 

The degree of the best model is: 2 

The number of models that were reduced is 76 

The total time this took was 667.580 

