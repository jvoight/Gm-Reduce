SetClassGroupBounds("GRH");

IdealShortVectorsProcess:=function(I, l, u : Minkowski:=true, timeout:=2);
//l,u are size in which to search over lattice, scaled by a medium sized vector in parallelepiped.
  if Degree(NumberField(Order(I))) gt 1 then
    if Minkowski eq true then
      Ibasis:=Basis(I);
      IxDen:=&*[ pp[1]^pp[2] : pp in Factorization(I*Denominator(I)) ];
      IxDen_gens,mxDen:=MinkowskiLattice(IxDen);
      Igens:=IxDen_gens/Denominator(I);
      BI:=Basis(Igens);
      dimL:=#BI;
      prec:=100+6*dimL;
      Rprec:=RealField(prec);
      VS:=VectorSpace(Rprec,dimL);
      BI_gens:=[ [ Rprec!v : v in Eltseq(gen) ] : gen in BI ];

      while Dimension(sub< VS | BI_gens >) lt dimL do
        prec:=prec+20;
        Rprec:=RealField(prec);
        VS:=VectorSpace(Rprec,dimL);
        BI_gens:=[ [ Rprec!v : v in Eltseq(gen) ] : gen in BI ];
      end while;

      LWB:=LatticeWithBasis(RMatrixSpace(Rprec,dimL,dimL)!&cat(BI_gens)); //lower the precision for efficiency
      OLWB:=LWB; // save original lattice in case wanted
      LWB:=CoordinateLattice(LWB);
      min_vec:=Min(LWB);
      avg_vec:=Determinant(LWB)^(1/(dimL));
      //EC:=EnumerationCost(LWB20, 10*Sqrt(avg_vec));
      SVP:=ShortVectorsProcess(LWB, avg_vec*l, avg_vec*u);

      SV:=[];
      while not(IsEmpty(SVP)) do
        Append(~SV, NextVector(SVP));
      end while;

      t:=Realtime(); //set a timeout
      SVcoord :=[];
      for w in SV do
        while Realtime(t) lt timeout do
          Append(~SVcoord, [ Round(c) : c in Eltseq(w) ]);
        end while;
      end for;

      SIelts := [ &+[w[i]*Ibasis[i] : i in [1..#Ibasis]] : w in SVcoord ] cat [1];
      //assert something to make sure there was no precision error
      return SIelts;

    else
      Igens := LLLBasis(I);
      // assert [ A/Denominator(I) : A in LLLBasis(I*Denominator(I)) ] eq Igens;
      Zn :=StandardLattice(#Igens);
      SVP:=ShortVectorsProcess(Zn, Ceiling(l),Ceiling(u));
      SV:=[];
      while not(IsEmpty(SVP)) do
        Append(~SV, NextVector(SVP));
      end while;
      MM := [ Eltseq(w) : w in SV ];
      SIelts:=[ (&+[ M[i]*Igens[i] : i in [1..#Igens] ]) : M in MM ] cat [1];
      return SIelts;
    end if;

  else
    tr, I_pr:=IsPrincipal(I);
    return [I_pr];
  end if;

end function;

function SmallFunctions(Qs, d)

  //        Qs, points which are "small"
  //        d, a degree bound; what d does genus 2 need?
  // Output: functions supported on {Qs} of degree <= d

  //The "small" points are the support of phi, phi-1.
  //sort "small" points by degree
  Qs := Sort(Qs, func<Q,R | Degree(Q)-Degree(R)>);

  Ds := [];
  QDs := [Divisor(Q) : Q in Qs]; //Make the points into divisors
  newDs := [Universe(QDs) | 0];

  //add divisors so that degree is <= d and collect in Ds
  while #newDs ne 0 do
    frontierDs := [];
    for newD in newDs do
      dnewD := Degree(newD);
      for RD in QDs do
        if Degree(RD) + dnewD le d then
          Append(~Ds, RD + newD);
          Append(~frontierDs, RD + newD);
        end if;
      end for;
    end for;
    newDs := frontierDs;
  end while;
  Ds:=Setseq(Set(Ds));

  xs := [];
  for Dden in Ds do
    for Dnum in Ds do

      D := Dden-Dnum;
      if D eq Parent(D)!0 then continue; end if;
      RR, mRR := RiemannRochSpace(D);
      if Dimension(RR) eq 1 then
        x := mRR(RR.1);
        if Divisor(x) ne Parent(Divisor(x))!0 then
          Append(~xs, x);
        end if;
      end if;
    end for;
  end for;
  //looks like every x appears twice
  return Setseq(Set(xs));
  //functions supported on the points of small degree as
  //genereated by the Riemann Roch space
end function;
//What is the difference in size from a divisor D and -D?
//Whats the relationship between the size of the places, the multiplicity and the size of f?

model := function(phi,x_op);
  //add 1/phi etc in here.
  fu := MinimalPolynomial(phi);
  fv := MinimalPolynomial(x_op);
  K := BaseRing(BaseRing(Parent(phi)));
  Kuvz<vz,uz,z> := PolynomialRing(K, 3);
  _<vf,uf,zf> := FieldOfFractions(Kuvz);

  fues := Eltseq(fu);
  fves := Eltseq(fv);
  fu := Numerator(&+[Evaluate(fues[i],zf)*uf^(i-1) : i in [1..#fues]]);
  fv := Numerator(&+[Evaluate(fves[i],zf)*vf^(i-1) : i in [1..#fves]]);

  fuv := Resultant(fu,fv,z);
  //groebner basis? 1/phi here etc
  /*
  _<u> := PolynomialRing(K);
  _<v> := PolynomialRing(Parent(u));
  cuv := Coefficients(fuv);
  muv := Monomials(fuv);
  return &+[cuv[i]*Evaluate(muv[i],[v,u,0]) : i in [1..#cuv]];
  */
  fuvFact := Factorization(fuv);
  if #fuvFact gt 1 then
    for j := 1 to #fuvFact do
      if Evaluate(fuvFact[j][1], [x_op, phi,0]) eq 0 then
        fuv := fuvFact[j][1];
        assert &and[Evaluate(fuvFact[k][1], [x_op, phi,0]) ne 0 : k in [j+1..#fuvFact]];
      end if;
    end for;
  end if;
  //_<u,v> := PolynomialRing(K,2);
  //return Evaluate(fuv,[v,u,0]);
  _<t,x> := PolynomialRing(K,2);
  return Evaluate(fuv,[x,t,0]);
  //x=v, t=u
end function;

PlaneModel := function(phi,x_op);
  return model(phi,x_op);
end function;

function PlaneModelGroebner(phi, x_op)
  //{Given a Belyi map phi, return a plane model for its domain such that t is the Belyi map}
  KC := Parent(phi);
  C := Curve(KC);
  K := BaseRing(C);
  nu := K.1;
  phi0 := Numerator(phi);
  phioo := Denominator(phi);
  x_op0 := Numerator(x_op);
  x_opoo := Denominator(x_op);
  if Genus(C) eq 0 then // defined on PP^1, so no curve equation
    R<X,Phi,X_op,u,v> := PolynomialRing(K,5);
    h := hom< KC -> R | [X]>;
    I := ideal< R | h(phioo)*Phi - h(phi0), h(x_opoo)*X_op - h(x_op0), v*h(phioo) - 1, u*h(x_opoo) - 1>; // need last equations to avoid points where phioo = 0
    // eliminate v to obtain plane equation
    basis := Basis(EliminationIdeal(I,{X_op,Phi}));
    assert #basis eq 1;
    new_eqn := basis[1];
    S<x,t> := PolynomialRing(K,2);
    h_plane := hom< Parent(new_eqn) -> S | [x,t,0,0,0] >;
  else // so the curve actually has an equation
    KC<x,y> := KC;
    R<X,Y,Phi,X_op,u,v> := PolynomialRing(K,6);
    h := hom< KC -> R | [X,Y]>;
    curve_eqn := DefiningEquation(AffinePatch(C,1));
    h_curve := hom< Parent(curve_eqn) -> R | [X,Y] >;
    // need last equation to avoid points where phioo = 0
    I := ideal< R | h_curve(curve_eqn), h(phioo)*Phi - h(phi0), h(x_opoo)*X_op - h(x_op0), v*h(phioo) - 1, u*h(x_opoo) - 1>; // need last equations to avoid points where phioo = 0
    // eliminate X and v to obtain plane equation
    basis := Basis(EliminationIdeal(I,{X_op,Phi}));
    assert #basis eq 1;
    new_eqn := basis[1];
    printf "new equation = %o\n", new_eqn;
    S<x,t> := PolynomialRing(K,2);
    h_plane := hom< Parent(new_eqn) -> S | [0,0,t,x,0,0] >;
    //h_plane := hom< Parent(new_eqn) -> S | X_op :-> x, Phi :-> t, X :-> 0, Y :-> 0, u :-> 0, v :-> 0 >;
  end if;
  return h_plane(new_eqn);
end function;

//What's the relationship between #model and #reduced_model? pretty much linear it would seem

//Pord<t1,x1>:=PolynomialRing(BaseRing(Parent(fuv)),2);
//Pord2<t2,x2>:=PolynomialRing(K,2: lex);

CoefficientValuationsSum:=function(f,pp)
  return &+[ Valuation(a,pp) : a in Coefficients(f) ],
  [ Valuation(a,pp) : a in Coefficients(f) ];
end function;

/*
  ZK:=Integers(K);
  LZK:=LLL(ZK);
  P5<b1,b2,b3,b4,b5>:=PolynomialRing(Rationals(),5);
  P5vars:=[b1,b2,b3,b4,b5];
  fuv_coefs:=Coefficients(fuv);
  [ &+[ P5vars[i]*cc[i] : i in [1..#P5vars] ] : cc in [ LZK!d : d in fuv_coefs ] ];
  P5xt<x,t>:=PolynomialRing(P5,2);
*/

BelyiObjectiveFunction:=function(fuv)
  K := BaseRing(Parent(fuv));
  u := Parent(fuv).1;
  v := Parent(fuv).2;
  Rx3<x1,x2,x3>:=PolynomialRing(Rationals(),3);
  mexps := [ Exponents(m) : m in Monomials(fuv) ];
  coefs:=Coefficients(fuv);
  assert &+[ coefs[i]*(u^mexps[i,1])*v^mexps[i,2] : i in [1..#mexps] ] eq fuv;
  return (&+[ m[1] : m in mexps])*x1 + (&+[ m[2] : m in mexps])*x2 + #mexps*x3;
end function;

MultivariateToUnivariate:=function(f)
  //turns an element in K[x,t] into an element K[x][t]
  fstring:=Sprint(f);
  K<nu1> := BaseRing(Parent(f));
  Kx<x>:=PolynomialRing(K);
  Kxt<t>:=PolynomialRing(Kx);
  return eval(fstring);
end function;


MonicToIntegral:=function(f);
  //scale the monic univariate polynomial to be integral
  assert IsMonic(f);
  K:=BaseRing(Parent(f));
  ZK:=Integers(K);
  coefs:=[ a : a in Coefficients(f) | a ne 0 ];
  pps:=[ a[1] : a in Factorization(&*coefs*ZK) ];
  scale:=1*ZK;
  for pp in pps do
    vals:=[ Valuation(cc,pp) : cc in coefs ];
    min:=Min(vals);
    scale:=scale*pp^(-min);
  end for;

  aprin, a := IsPrincipal(scale);
  if aprin eq false then a:=IdealShortVectorsProcess(scale, 0.00001, 2: Minkowski:=false)[1]; end if;
  f_new:=f*a;

  return f_new, 1/a;
end function;



PolynomialToFactoredString:=function(f)
  //factorise the polynomial and return it as a string
  //Needs to be a multivariate polynomial in K[x][t]

  coefs:=Coefficients(f);
  mon:=Monomials(f);
  str:="";
  for j in [1..#coefs] do

    if j ne 1 then
      str:=str cat " + ";
    end if;

    if Degree(coefs[j]) ne 0 then
      a:=LeadingCoefficient(coefs[j]);
      fac:=Factorization(coefs[j]);
      list:=[];
      for item in fac do
        int,co:= MonicToIntegral(item[1]);
        Append(~list,<co,int,item[2]>);
      end for;

      c:=a*(&*[b[1] : b in list]);
      str:= str cat Sprintf("(%o)",c);
      for i in [1..#list] do
        str:= str cat Sprintf("*(%o)^%o", list[i,2],list[i,3]);
      end for;

    else
      str:=str cat Sprint(coefs[j]);
    end if;

    if j ne 1 then
      str:=str cat Sprintf("*%o",mon[j]);
    end if;

  end for;

  K<nu1>:=BaseRing(BaseRing(f));
  Kx<x>:=PolynomialRing(K);
  Kxt<t>:=PolynomialRing(Kx);

  assert f eq eval(str);
  return str;
end function;



reducemodel_padic := function(f : Polyhedron:=false);
  //input: a multivariate polynomial f \in K[z_1,..,z_n]
  //output: minimal and integral c*f(a_1z_1,...,a_nz_n) and [a_1,...,a_n,c]
  K := BaseRing(Parent(f));
  variables:=[ Parent(f).i : i in [1..#Names(Parent(f))] ];
  n:=#variables;
  ZK := Integers(K);
  k:=Integers();

  //Rx3<x1,x2,x3>:=PolynomialRing(Rationals(),3);
  //obj_fun:= BelyiObjectiveFunction(fuv);

  coefs_and_monomials:= [ [Coefficients(f)[i],Monomials(f)[i]] : i in [1..#Coefficients(f)] | Coefficients(f)[i] ne 0 ];
  mexps := [ Exponents(m[2]) : m in coefs_and_monomials ];
  m:=#mexps;
  coefs:=[ K!a[1] : a  in coefs_and_monomials ];
  //assert &+[ coefs[i]*(u^mexps[i,1])*v^mexps[i,2] : i in [1..#mexps] ] eq fuv;
  obj_coefs:= [ &+[ m[i] : m in mexps] : i in [1..n] ] cat [m];
  obj := Matrix(k,1,n+1, obj_coefs);

  //S is the prime divisors of all norms of numerators and denominators of coeffients
  S := &cat[TrialDivision(Integers()!Norm(Numerator(s))) : s in Coefficients(f) | s ne 0 ]
		 cat &cat[TrialDivision(Integers()!Norm(Denominator(s))) : s in Coefficients(f) | s ne 0 ];
  //do we need a bound for trial division?
  S := SequenceToSet([s[1] : s in S]);
  SS:=&cat[ [pp[1] : pp in Factorization(p*Integers(K))] : p in S ];

  rescaling_ideals:=[[ 1*ZK : i in [1..n+1] ]];
  new_c := 1;
  new_f:=f;

  for pp in SS do
	  cvals := [ Valuation(c,pp) : c in coefs  ];
    //valuations at this prime pp of the coefficients
    lhs_coefs:= [ [ mexps[i,j] : j in [1..n] ] cat [1] : i in [1..m] ];
    lhs := Matrix(k, lhs_coefs); //constraints
    rel := Matrix(k,[[1] : ef in mexps]);             //lhs greater than rhs
    rhs := Matrix(k, [[-cf] : cf in cvals]);          //valuations

    halfspaces:=[ HalfspaceToPolyhedron(Eltseq(Rows(lhs)[i]),Eltseq(rhs)[i]) : i in [1..#Rows(lhs)] ];
    poly:= &meet[ POL : POL in halfspaces ];
    //find the minimum of the objective function in the region, either using integral vertices or the linear program
    if Polyhedron eq false then
      L := LPProcess(k, n+1);
      SetObjectiveFunction(L, obj);
      AddConstraints(L, lhs, rhs : Rel := "ge");
      //UnsetBounds(L) doesn't work
      //These are lower bounds on the solution
      for i in [1..n+1] do  SetLowerBound(L, i, k!-10000); end for;
      soln,state:=Solution(L);
      //ProfilePrintByTotalTime(:Max:=40);
      assert state eq 0;
      min:=EvaluateAt(L,soln);

      //assert [Eltseq(r) : r in Rays(InfinitePart(int_poly_old))] eq [[0,1],[1,0]] or
      //&and[&and[IsIntegral(cvt) : cvt in Eltseq(vt)] : vt in vts_old];
      //vts := [ [ Ceiling(cvt) : cvt in Eltseq(vt)] : vt in vts];
      //ceiling?
    //else
      //int_poly := IntegralPart(poly);
      //vts:=Vertices(int_poly);

      //vertex_solns:=[];
      //for vt in vts do
        //Append(~vertex_solns, Evaluate(obj_fun,Eltseq(vt)));
      //end for;
      //min:=Minimum(vertex_solns);
    end if;

    //Now we intersect our polyhedron with the 'plane of minimal solutions'
    minimal_hyperplane:=HyperplaneToPolyhedron(Eltseq(obj),min);
    poly := poly meet minimal_hyperplane;
    int_poly := IntegralPart(poly);
    assert IsEmpty(poly) eq false;
    assert IsPolytope(poly);
    min_points:=Setseq(Points(int_poly));
    //all of the points at which the objective function is minimal.

    //assert Eltseq(soln) in [ Eltseq(a) : a in min_points ];

    //find the points which give something principal
    Cl,h:=ClassGroup(K); hin:=Inverse(h);
    prin:=Order(hin(pp));
    principal_points := [ vv : vv in min_points | [ IsDivisibleBy(a,prin) : a  in Eltseq(vv) ] eq [true,true,true] ];

    if principal_points eq [] then
      points_loop := min_points; principal:=false;
    else
      points_loop := principal_points; principal:=true;
    end if;

    //all triples of ideals to try rescaling by
    rescaling_ideals:=&cat[ [ [ (ideals[i])*(pp^Eltseq(pt)[i]) : i in [1..#Eltseq(pt)] ] : pt in points_loop ] : ideals in rescaling_ideals ];

  end for;

  new_fuvs:=[];
  for vv in rescaling_ideals do
    scaling_factors:= [ ];
    for w in vv do
      aprin, a := IsPrincipal(w);
      if aprin eq false then a:=IdealShortVectorsProcess(w, 0.00001, 2: Minkowski:=false)[1]; end if;
      Append(~scaling_factors, a);
    end for;

    if n eq 1 then
      guv:=Evaluate(new_f,scaling_factors[1]*variables[1])*scaling_factors[n+1];
    else
      guv:=Evaluate(new_f,[scaling_factors[i]*variables[i] : i in [1..n]])*scaling_factors[n+1];
    end if;
    Append(~new_fuvs, <#Sprint(guv),guv,scaling_factors>);
  end for;

  Sort(~new_fuvs);
  new_fuv:=new_fuvs[1,2];
  new_scaling:= new_fuvs[1,3];

  return new_fuv, new_scaling;
end function;

reducemodel_units := function(fuv : Polyhedron:=false);
  K := BaseRing(Parent(fuv));
  u := Parent(fuv).1;
  v := Parent(fuv).2;
  ZK := Integers(K);
  k:=Integers();
  Rx3<x1,x2,x3>:=PolynomialRing(Rationals(),3);

  mexps := [ Exponents(m) : m in Monomials(fuv) ];
  coefs:=Coefficients(fuv);
  assert &+[ coefs[i]*(u^mexps[i,1])*v^mexps[i,2] : i in [1..#mexps] ] eq fuv;

  k1:=RealField(20);
  UK,mUK:=UnitGroup(K);
  M,phi:=MinkowskiSpace(K);
  UU:= [ K!(mUK(eps)) : eps in Generators(UK) | not(IsFinite(eps)) ];

    N:=#mexps*Dimension(M);
    L1 := LPProcess(k1, 3*#UU+N);
    obj1:=Matrix(k1,1,3*#UU+N,[0 : i in [1..3*#UU]] cat [1 : i in [1..N]]);
    SetObjectiveFunction(L1, obj1);

    for n in [1..#mexps] do
      alpha_norm:=Log(Abs(Norm(coefs[n])))/(Dimension(M));
      log_coef:= [ Log(Abs(alpha)) : alpha in Eltseq(phi(coefs[n])) ];
      for m in [1..Dimension(M)] do
        extra_var1:=[ 0 : k in [1..N-1] ];
        Insert(~extra_var1, (n-1)*Dimension(M) +m, -1);

        phi_eps:= [ Log(Abs(Eltseq(phi(e))[m])) : e in UU ];
        lhs1:= Matrix(k1,1,3*#UU+N, &cat[ [ mexps[n,1]*log_eps, mexps[n,2]*log_eps, log_eps] : log_eps in phi_eps ] cat extra_var1);
        lhs2:= Matrix(k1,1,3*#UU+N, &cat[ [ -mexps[n,1]*log_eps, -mexps[n,2]*log_eps, -log_eps] : log_eps in phi_eps ] cat extra_var1);

        rhs1:= Matrix(k1,[[-log_coef[m] + alpha_norm]]);
        rhs2:= Matrix(k1,[[log_coef[m] - alpha_norm]]);

        AddConstraints(L1, lhs1, rhs1 : Rel := "le");
        AddConstraints(L1, lhs2, rhs2 : Rel := "le");

        AddConstraints(L1, Matrix(k1,1,3*#UU+N, [0 : i in [1..3*#UU]] cat extra_var1), Matrix(k1,[[0]]) : Rel :="le");

      end for;
    end for;

    for i in [1..3*#UU+N] do
      SetLowerBound(L1, i, k1!-10000);
    end for;
    soln,state:=Solution(L1);
    //ProfilePrintByTotalTime(:Max:=40);
    assert state eq 0;

    eps_soln:=[ Eltseq(soln)[i] : i in [1..3*#UU] ];
    eps_rounded:= [ Integers()!Round(e) : e in eps_soln ];
    a:= &*[ UU[j]^eps_rounded[1+3*(j-1)] : j in [1..#UU] ];
    b:= &*[ UU[j]^eps_rounded[2+3*(j-1)] : j in [1..#UU] ];
    c:= &*[ UU[j]^eps_rounded[3+3*(j-1)] : j in [1..#UU] ];

    guv:=Evaluate(fuv,[a*u,b*v])*c;

    return guv, [a,b,c];
end function;

reducemodel_old := function(fuv : Polyhedron:=false);
  K := BaseRing(Parent(fuv));
  u := Parent(fuv).1;
  v := Parent(fuv).2;
  ZK := Integers(K);
  mexps := [ Exponents(m) : m in Monomials(fuv) ];
  coefs:=Coefficients(fuv);

  //S is the prime divisors of all norms of numerators and denominators of coeffients
  S := &cat[TrialDivision(Integers()!Norm(Numerator(s))) : s in Coefficients(fuv)]
		 cat &cat[TrialDivision(Integers()!Norm(Denominator(s))) : s in Coefficients(fuv)];
  //do we need a bound for trial division?
  S := SequenceToSet([s[1] : s in S]);
  aa := 1*ZK;
  bb := 1*ZK;
  cc := 1*ZK;
  for p in S do
	  for pp in [pp[1] : pp in Factorization(p*Integers(K))] do
	    cvals := [ Valuation(c,pp) : c in coefs ];
      //valuations at this prime pp of the coefficients
	    k := Integers();
	    nodemins := [];
	    for jc := 1 to #mexps do
        //run over the exponents of the coefficients
    		e0f0 := mexps[jc];
	     	e0,f0 := Explode(e0f0);
		    c0 := cvals[jc];
		    lhs := Matrix(k,[[ef[1]-e0, ef[2]-f0] : ef in mexps]); //The difference in the exponents of each coeffient and the exponents of this coefficient
		    rel := Matrix(k,[[1] : ef in mexps]);
		    rhs := Matrix(k,[[c0-c] : c in cvals]); //The difference between valuations
		    obj := Matrix(k,1,2,[&+[mexps[j][i]-e0f0[i] : j in [1..#mexps]] : i in [1..2]]);
		    // integer linear program doesn't work?

        if Polyhedron eq false then
	        soln:=MinimalIntegerSolution(lhs, rel, rhs, obj);
          soln:= Matrix(k,1,2,[ Round(soln[1,1]), Round(soln[1,2]) ]);
          objvals:=obj[1,1]*soln[1,1] + obj[1,2]*soln[1,2];
          Append(~nodemins, <objvals, soln, c0, e0f0>);
        else
          poly := PolyhedronWithInequalities(
  						  [ Eltseq(r) : r in Rows(lhs) ],
  						   ChangeUniverse(Eltseq(rhs),Integers()));
  		    // it's a union of cones and obj is linear, so min is at a vertex
  		    poly := IntegralPart(poly);
  		    vts := Vertices(poly);
  		    assert [Eltseq(r) : r in Rays(InfinitePart(poly))] eq [[0,1],[1,0]] or
  				    &and[&and[IsIntegral(cvt) : cvt in Eltseq(vt)] : vt in vts];
  		    vts := [ [ Ceiling(cvt) : cvt in Eltseq(vt)] : vt in vts];
          //looks like we might get a couple of vertices using this method
          //and only one with the linear program, which is better?
  		    if #vts gt 0 then
  		      objvals := [(obj*Matrix(k,2,1,Eltseq(vt)))[1,1] - #mexps*c0 : vt in vts];
  		      mval, mind := Min(objvals);
  		      Append(~nodemins, <mval, vts[mind], c0, e0f0>);
  		    end if;
        end if;
      end for;

      nodemins := Sort(nodemins, func<x,y | x[1]-y[1]>);
    	ef := Eltseq(nodemins[1][2]);
    	e0f0 := Eltseq(nodemins[1][4]);
    	aa *:= pp^(ef[1]);
    	bb *:= pp^(ef[2]);
    	cc *:= pp^(nodemins[1][3]+ef[1]*e0f0[1]+ef[2]*e0f0[2]);

    end for;
  end for;

  aprin, a := IsPrincipal(aa);
  if aprin eq false then a:=IdealShortVectorsProcess(aa, 0.00001, 2: Minkowski:=false)[1]; end if;
  bprin, b := IsPrincipal(bb);
  if bprin eq false then b:=IdealShortVectorsProcess(bb, 0.00001, 2: Minkowski:=false)[1]; end if;
  cprin, c := IsPrincipal(cc);
  if cprin eq false then c:=IdealShortVectorsProcess(cc, 0.00001, 2: Minkowski:=false)[1]; end if;
  //include units
/*
  cvals := [ Valuation(c,pp) : c in coefs ];
  //valuations at this prime pp of the coefficients
  k := Integers();
  nodemins := [];
  for jc := 1 to #mexps do
    //run over the exponents of the coefficients
    e0f0 := mexps[jc];
    e0,f0 := Explode(e0f0);
    c0 := cvals[jc];
    lhs := Matrix(k,[[ef[1]-e0, ef[2]-f0] : ef in mexps]); //The difference in the exponents of each coeffient and the exponents of this coefficient
    rel := Matrix(k,[[1] : ef in mexps]);
    rhs := Matrix(k,[[c0-c] : c in cvals]); //The difference between valuations
    obj := Matrix(k,1,2,[&+[mexps[j][i]-e0f0[i] : j in [1..#mexps]] : i in [1..2]]);
    // integer linear program doesn't work?

    if Polyhedron eq false then
      soln:=MinimalIntegerSolution(lhs, rel, rhs, obj);
      soln:= Matrix(k,1,2,[ Round(soln[1,1]), Round(soln[1,2]) ]);
      objvals:=obj[1,1]*soln[1,1] + obj[1,2]*soln[1,2];
      Append(~nodemins, <objvals, soln, c0, e0f0>);
*/

/*
  fuv_new:=Evaluate(fuv,[a*u,b*v])/c;
  UK,mUK:=UnitGroup(K);
  eps:=K!(mUK(UK.2));

  ff_units:=[];
  for i in [-10..10] do
    for j in [-10..10] do
       fU:=Evaluate(fuv_new,[u,v*eps^i])/eps^j;
       Append(~ff_units,<#Sprint(fU),fU>);
    end for;
  end for;
  ff_units:=Sort(ff_units);
*/
  return [ Evaluate(fuv,[a*u,b*v])/c, 1/c];
end function;
