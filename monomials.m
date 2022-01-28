AttachSpec("../endomorphisms/endomorphisms/magma/spec"); // have to change if endomorphisms repo is elsewhere
AttachSpec("../Belyi/Code/spec"); // have to change if Belyi repo is elsewhere
// includes intrinsic S3Action(tau, phi)

/*
  OK := Integers(K);
  dens := [];
  for el in Coefficients(Numerator(phi)) do
    Append(~dens, Denominator(el));
  end for;
  for el in Coefficients(Denominator(phi)) do
    Append(~dens, Denominator(el));
  end for;
  dens;
  cop := &*dens;
  cop *:= Discriminant(OK);
  primes := PrimesUpTo(100, K : coprime_to := cop);
*/

intrinsic ReduceRationalFunction(X::Crv, phi::FldFunFracSchElt, P::RngOrdIdl) -> Any
  {Given a Belyi map phi defined on a curve X and a prime ideal of the field of definition of X, return the reductions of X and phi mod P}

  // make curve over finite field
  KX := Parent(phi);
  K := BaseRing(X);
  OK := Integers(K);
  FF, res_mp := ResidueClassField(P);
  X_FF := Reduction(X, P);
  X_FF := Curve(X_FF);
  KX_FF := FunctionField(X_FF);

  AX := CoordinateRing(X);
  N_gens := #GeneratorsSequence(AX) - 1; // number of generators for the function field
  KX_gens := [];
  KX_FF_gens := [];
  for i := 1 to N_gens do
    Append(~KX_gens, KX.i);
    Append(~KX_FF_gens, KX_FF.i);
  end for;

  // reduce map 
  num_cs, num_mons := CoefficientsAndMonomials(Numerator(phi));
  num_pows := [Exponents(el) : el in num_mons];
  den_cs, den_mons := CoefficientsAndMonomials(Denominator(phi));
  den_pows := [Exponents(el) : el in den_mons];
  phi_FF_num := KX_FF!0;
  for i := 1 to #num_cs do
    phi_FF_num +:= res_mp(num_cs[i])*&*[KX_FF.j^num_pows[i][j] : j in [1..N_gens]];
  end for;
  phi_FF_den := KX_FF!0;
  for i := 1 to #den_cs do
    phi_FF_den +:= res_mp(den_cs[i])*&*[KX_FF.j^den_pows[i][j] : j in [1..N_gens]];
  end for;
  return X_FF, phi_FF_num/phi_FF_den;
end intrinsic;

//SortSmallFunctions
